"""Shared calibrated reader for Fig. 4D (rescaled bundle contraction master collapse)."""

from __future__ import annotations

import numpy as np
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    xt = C.match_ticks_regular(xt, 4, tol=4.0)
    yt = C.find_ticks(dark, vr0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=-1, band=(8, 45))
    yt = C.match_ticks_regular(yt, 6, tol=4.0)

    xa = C.verify_axis("linear", [(vl0[0], 0.0)] + list(zip(xt, [1.0, 2.0, 3.0, 4.0])) + [(vr0[0], 5.0)], 1.5)
    ya = C.verify_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), 1.5)
    return xa, ya, int(ht[2]) + 10, int(hb[1]) - 10, int(vl0[2]) + 10, int(vr0[1]) - 10


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)

    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    # Exclude legend box in upper-right (t/tau > 3.1, y > 0.60)
    inb[: int(ya.to_pixel(0.60)), int(xa.to_pixel(3.1)) :] = False

    if which == "2p5pct":
        mask = inb & (g - r > 45.0) & (g - b > 45.0) & (g > 110.0)
        label = "2.5% green"
    elif which == "5pct":
        mask = inb & (r - g > 50.0) & (r - b > 45.0) & (r > 160.0)
        label = "5.0% red"
    elif which == "10pct":
        mask = inb & (b - r > 50.0) & (b > 140.0) & (r < 160.0)
        label = "10.0% blue"
    else:
        raise ValueError(which)

    W = mask.shape[1]
    out = []
    for x in at:
        xf = float(x)
        cf = float(xa.to_pixel(xf))
        half_dx = 0.12 if xf >= 0.4 else 0.07
        c_lo = max(c0, int(round(xa.to_pixel(xf - half_dx))))
        c_hi = min(c1, int(round(xa.to_pixel(xf + half_dx))))
        rr, cc = np.nonzero(mask[:, c_lo : c_hi + 1])
        if len(rr) < 15:
            out.append({"x": xf, "y": None, "found": False, "note": f"no {label} points near t/tau={xf}"})
            continue
        xs_loc = cc.astype(float) + c_lo + 0.5 - cf
        ys_loc = rr.astype(float) + 0.5
        if np.ptp(xs_loc) > 8.0:
            slope, r_est = np.polyfit(xs_loc, ys_loc, 1)
            resid = ys_loc - (slope * xs_loc + r_est)
            r_p10, r_p90 = float(r_est + np.percentile(resid, 10)), float(r_est + np.percentile(resid, 90))
        else:
            r_est = float(np.median(ys_loc))
            r_p10, r_p90 = float(np.percentile(ys_loc, 10)), float(np.percentile(ys_loc, 90))
        pad = 2.0 + ya.residual_px
        yv = float(ya.to_value(r_est))
        ya_, yb_ = float(ya.to_value(min(r_p10, r_est) - pad)), float(ya.to_value(max(r_p90, r_est) + pad))
        out.append(
            {
                "x": xf,
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "px": [cf, float(r_est)],
                "px_lo": float(min(r_p10, r_est) - pad),
                "px_hi": float(max(r_p90, r_est) + pad),
                "note": f"Fig4D {label} cloud 10-90% (n={len(rr)}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            }
        )
    return out
