"""Shared calibrated reader for Fig. 2A (tactoid aspect ratio L/r and anchoring strength omega vs filamin %)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, thresh=150.0, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    xt = C.match_ticks_regular(xt, 6, tol=4.0)
    yt_l = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(8, 45))
    yt_l = C.match_ticks_regular(yt_l, 5, tol=4.0)
    yt_r = C.find_ticks(dark, vr0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=-1, band=(8, 45))
    yt_r = C.match_ticks_regular(yt_r, 5, tol=4.0)

    xa = C.verify_axis("linear", [(vl0[0], 2.0)] + list(zip(xt, [4.0, 6.0, 8.0, 10.0, 12.0, 14.0])) + [(vr0[0], 16.0)], 1.0)
    ya_l = C.verify_axis("linear", [(ht[0], 3.5)] + list(zip(yt_l, [3.0, 2.5, 2.0, 1.5, 1.0])) + [(hb[0], 0.5)], 1.0)
    ya_r = C.verify_axis("linear", [(ht[0], 3.5)] + list(zip(yt_r, [3.0, 2.5, 2.0, 1.5, 1.0])) + [(hb[0], 0.5)], 1.0)
    return xa, ya_l, ya_r, int(ht[2]) + 10, int(hb[1]) - 10, int(vl0[2]) + 40, int(vr0[1]) - 40


def read_series(image: np.ndarray, which: str) -> list[dict]:
    xa, ya_l, ya_r, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True

    if which == "omega":
        green = inb & (g - r > 40) & (g > 100)
        # Fast distance transform locates solid diamond cores (half-diagonal ~40 px) vs thin error-bar stems (~4 px)
        cores = ndimage.distance_transform_edt(green) >= 15.0
        lab, n = ndimage.label(cores)
        out = []
        for i in range(1, n + 1):
            rr, cc = np.nonzero(lab == i)
            if len(rr) < 200:
                continue
            c_m, r_m = float(cc.mean()) + 0.5, float(rr.mean()) + 0.5
            bar_lo, bar_hi, bar_lo_out, bar_hi_out, row_hi, row_lo = C.trace_dashed_error_bar(
                green, c_m, int(rr.min()), int(rr.max()), r0, r1, ya_r
            )
            half_px = 3.0 + ya_r.residual_px
            yv = float(ya_r.to_value(r_m))
            ya_, yb_ = float(ya_r.to_value(r_m - half_px)), float(ya_r.to_value(r_m + half_px))
            out.append(
                {
                    "x": float(xa.to_value(c_m)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "bar_lo": bar_lo,
                    "bar_hi": bar_hi,
                    "bar_lo_outer": bar_lo_out,
                    "bar_hi_outer": bar_hi_out,
                    "px": [c_m, r_m],
                    "px_lo": float(row_hi),
                    "px_hi": float(row_lo),
                    "note": f"Fig2A green diamond, bar=[{bar_lo:.3f},{bar_hi:.3f}], cal_res=({xa.residual_px:.2f},{ya_r.residual_px:.2f})px",
                }
            )
        out.sort(key=lambda d: d["x"])
        return out

    if which == "aspect_ratio":
        blk = inb & (r < 60) & (g < 60) & (b < 60)
        # Exclude inset tactoid schematic (x < 7.2%, y < 1.45)
        blk[int(ya_l.to_pixel(1.45)) :, : int(xa.to_pixel(7.2))] = False
        holes = ndimage.binary_fill_holes(blk) & ~blk
        lab, n = ndimage.label(holes)
        out = []
        for i in range(1, n + 1):
            rr, cc = np.nonzero(lab == i)
            if len(rr) < 1500:
                continue
            c_m, r_m = float(cc.mean()) + 0.5, float(rr.mean()) + 0.5
            bar_lo, bar_hi, bar_lo_out, bar_hi_out, row_hi, row_lo = C.trace_dashed_error_bar(
                blk, c_m, int(rr.min()) - 20, int(rr.max()) + 20, r0, r1, ya_l
            )
            half_px = 3.0 + ya_l.residual_px
            yv = float(ya_l.to_value(r_m))
            ya_, yb_ = float(ya_l.to_value(r_m - half_px)), float(ya_l.to_value(r_m + half_px))
            out.append(
                {
                    "x": float(xa.to_value(c_m)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "bar_lo": bar_lo,
                    "bar_hi": bar_hi,
                    "bar_lo_outer": bar_lo_out,
                    "bar_hi_outer": bar_hi_out,
                    "px": [c_m, r_m],
                    "px_lo": float(row_hi),
                    "px_hi": float(row_lo),
                    "note": f"Fig2A open circle, bar=[{bar_lo:.3f},{bar_hi:.3f}], cal_res=({xa.residual_px:.2f},{ya_l.residual_px:.2f})px",
                }
            )
        out.sort(key=lambda d: d["x"])
        return out

    raise ValueError(which)
