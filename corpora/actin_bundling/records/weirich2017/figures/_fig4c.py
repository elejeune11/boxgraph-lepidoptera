"""Shared calibrated reader for Fig. 4C (bundle contraction L(t)/L_i vs time for 2.5% and 10% filamin)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    xt = C.match_ticks_regular(xt, 2, tol=5.0)
    yt = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(8, 45))
    yt = C.match_ticks_regular(yt, 5, tol=4.0)

    # Interior x-ticks 500, 1000 plus left spine 0 give exact <0.2 px calibration (right spine is +4 px past 1500)
    xa = C.verify_axis("linear", [(vl0[0], 0.0), (xt[0], 500.0), (xt[1], 1000.0)], max_residual_px=1.0)
    ya = C.verify_axis("linear", list(zip(yt, [0.9, 0.8, 0.7, 0.6, 0.5])), max_residual_px=1.0)
    return xa, ya, int(ht[2]) + 10, int(hb[1]) - 10, int(vl0[2]) + 10, int(vr0[1]) - 10


def _sample_marker_trace(mask: np.ndarray, xa, ya, at: list[float], half_dt: float, label: str) -> list[dict]:
    """Sample a densely plotted single-trace experimental marker series over narrow +-half_dt bins."""
    H, W = mask.shape
    out = []
    for x in at:
        xf = float(x)
        cf = float(xa.to_pixel(xf))
        c_lo = max(0, int(round(xa.to_pixel(xf - half_dt))))
        c_hi = min(W - 1, int(round(xa.to_pixel(xf + half_dt))))
        sub = mask[:, c_lo : c_hi + 1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 20:
            out.append({"x": xf, "y": None, "found": False, "note": f"no {label} pixels near t={xf}"})
            continue
        r_med = float(np.median(rr)) + 0.5
        # Use local linear fit if window spans multiple markers so slope across bin does not bias center
        xs_loc = cc.astype(float) + c_lo + 0.5 - cf
        ys_loc = rr.astype(float) + 0.5
        if np.ptp(xs_loc) > 10.0:
            slope, intercept = np.polyfit(xs_loc, ys_loc, 1)
            r_est = float(intercept)
            resid = np.abs(ys_loc - (slope * xs_loc + intercept))
            spread = float(np.percentile(resid, 75))
        else:
            r_est = r_med
            slope = 0.0
            spread = float(np.percentile(np.abs(ys_loc - r_med), 75))
        half_px = max(spread, 6.0) + ya.residual_px + abs(slope) * (1.0 + xa.residual_px)
        yv = float(ya.to_value(r_est))
        ya_, yb_ = float(ya.to_value(r_est - half_px)), float(ya.to_value(r_est + half_px))
        out.append(
            {
                "x": xf,
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "px": [cf, r_est],
                "px_lo": float(r_est - half_px),
                "px_hi": float(r_est + half_px),
                "note": f"Fig4C {label} (n={len(rr)}px, hw={half_px:.1f}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            }
        )
    return out


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)

    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    # Exclude legend box in top-right (t > 1080, y > 0.84)
    inb[: int(ya.to_pixel(0.84)), int(xa.to_pixel(1080.0)) :] = False

    if which == "2p5pct_data":
        # Red open circles
        red = inb & (r - g > 80.0) & (r - b > 80.0)
        return _sample_marker_trace(red, xa, ya, at, half_dt=22.0, label="2.5% red circles")

    if which == "10pct_data":
        # Blue filled circles
        blue = inb & (b - r > 70.0) & (b > 140.0)
        return _sample_marker_trace(blue, xa, ya, at, half_dt=18.0, label="10.0% blue circles")

    if which in ("2p5pct_fit", "10pct_fit"):
        # Pure neutral black dashes; inset 78 px from all spines to completely exclude inward ticks
        blk = (r < 45) & (g < 45) & (b < 45) & (np.abs(r - g) < 15) & (np.abs(r - b) < 15)
        m = np.zeros_like(blk)
        r_lo, r_hi = r0 + 75, r1 - 75
        c_lo, c_hi = c0 + 75, c1 - 75
        m[r_lo:r_hi, c_lo:c_hi] = blk[r_lo:r_hi, c_lo:c_hi]
        m[: int(ya.to_pixel(0.84)), int(xa.to_pixel(1080.0)) :] = False
        # Mask L_i schematic (t < 240, y < 0.76) and L_f schematic (t in 1010..1240, y < 0.52)
        m[int(ya.to_pixel(0.76)) :, : int(xa.to_pixel(240.0))] = False
        m[int(ya.to_pixel(0.52)) :, int(xa.to_pixel(1010.0)) : int(xa.to_pixel(1240.0))] = False

        # Separate upper (2.5%) vs lower (10.0%) dashes cleanly using smooth midpoint separator y_split(t)
        lab, _ = ndimage.label(m)
        curve_mask = np.zeros_like(m)
        for bl in C.find_blobs(m, min_area=55, max_area=2000):
            if bl["h"] > 90 or bl["w"] > 120:
                continue
            tv = float(xa.to_value(bl["c"]))
            yv = float(ya.to_value(bl["r"]))
            y_split = 0.635 + 0.365 * np.exp(-tv / 230.0)
            if (which == "2p5pct_fit" and yv > y_split) or (which == "10pct_fit" and yv <= y_split):
                sl = (slice(bl["r0"], bl["r1"] + 1), slice(bl["c0"], bl["c1"] + 1))
                curve_mask[sl] |= lab[sl] == bl["label"]
        return [
            C.sample_dashed_curve_y(
                curve_mask, xa, ya, float(x), r_lo, r_hi, max_search_dx=220, sigma_dx=90.0, deg=2, note_prefix=f"Fig4C {which}: "
            )
            for x in at
        ]

    raise ValueError(which)
