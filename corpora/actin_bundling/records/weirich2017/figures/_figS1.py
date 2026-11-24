"""Shared calibrated reader for SI Fig. S1 (bleached & unbleached FRAP fluorescence intensity vs time)."""

from __future__ import annotations

import numpy as np
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45), max_width=30)
    xt = C.match_ticks_regular(xt, 5, tol=4.0)
    yt = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(8, 45), max_width=30)
    yt = C.match_ticks_regular(yt, 4, tol=4.0)

    xa = C.verify_axis("linear", list(zip(xt, [0.0, 200.0, 400.0, 600.0, 800.0])), max_residual_px=1.5)
    ya = C.verify_axis(
        "linear", [(ht[0], 1200.0)] + list(zip(yt, [1000.0, 800.0, 600.0, 400.0])) + [(hb[0], 200.0)], max_residual_px=1.5
    )
    return xa, ya, int(ht[2]) + 15, int(hb[1]) - 15, int(vl0[2]) + 15, int(vr0[1]) - 15


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)

    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    # Exclude bottom-right inset micrograph (t > 620 s, intensity < 540)
    inb[int(ya.to_pixel(540.0)) :, int(xa.to_pixel(620.0)) :] = False

    if which == "bleached":
        # Pink/magenta dots (bleached side, RGB ~(240, 85, 190))
        mask = inb & (r - g > 65.0) & (r > 180.0) & (b > 110.0)
        # Note: exclude the 2 pre-bleach pink dots at t <= 0, I > 1050 when sampling post-bleach t >= 0
        mask[: int(ya.to_pixel(1000.0)), :] = False
        label = "bleached (pink)"
    elif which == "unbleached":
        # Dark blue dots (unbleached side, RGB ~(10, 20, 200))
        mask = inb & (b - r > 80.0) & (b > 130.0)
        label = "unbleached (blue)"
    else:
        raise ValueError(which)

    out = []
    for x in at:
        xf = float(x)
        cf = float(xa.to_pixel(xf))
        half_dt = 16.0 if xf > 10.0 else 12.0
        t_lo = max(-2.0, xf - half_dt) if xf <= 5.0 else xf - half_dt
        c_lo = max(c0, int(round(xa.to_pixel(t_lo))))
        c_hi = min(c1, int(round(xa.to_pixel(xf + half_dt))))
        rr, cc = np.nonzero(mask[:, c_lo : c_hi + 1])
        if len(rr) < 20:
            out.append({"x": xf, "y": None, "found": False, "note": f"no {label} dots near t={xf} s"})
            continue
        xs_loc = cc.astype(float) + c_lo + 0.5 - cf
        ys_loc = rr.astype(float) + 0.5
        if np.ptp(xs_loc) > 10.0:
            slope, r_est = np.polyfit(xs_loc, ys_loc, 1)
            resid = ys_loc - (slope * xs_loc + r_est)
            r_lo_p = float(r_est + np.percentile(resid, 15))
            r_hi_p = float(r_est + np.percentile(resid, 85))
        else:
            r_est = float(np.median(ys_loc))
            r_lo_p, r_hi_p = float(np.percentile(ys_loc, 15)), float(np.percentile(ys_loc, 85))
        pad = 3.0 + ya.residual_px
        yv = float(ya.to_value(r_est))
        ya_, yb_ = float(ya.to_value(min(r_lo_p, r_est) - pad)), float(ya.to_value(max(r_hi_p, r_est) + pad))
        out.append(
            {
                "x": xf,
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "px": [cf, float(r_est)],
                "px_lo": float(min(r_lo_p, r_est) - pad),
                "px_hi": float(max(r_hi_p, r_est) + pad),
                "note": f"FigS1 {label} (n={len(rr)}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            }
        )
    return out
