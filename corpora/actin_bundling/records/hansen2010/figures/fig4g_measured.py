"""Calibrated reader for Fig. 4G (p. 7): barbed-end polymerization rate vs [VASP]tet with SD error bars."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig4g(image)
    cL, rT, cR, rB = fr["box"]
    blk = C.black_ink_mask(image, gray_max=110.0, sat_max=35.0)
    # Exclude Kd annotation box (x > 40, y < 12)
    blk[int(ya.to_pixel(12.0)) :, int(xa.to_pixel(40.0)) :] = False

    # Detect interior filled markers (1.25..100 nM) via morphological opening
    sub = np.zeros_like(blk)
    sub[int(rT) + 10 : int(rB) - 10, int(cL) + 4 : int(cR) - 10] = blk[int(rT) + 10 : int(rB) - 10, int(cL) + 4 : int(cR) - 10]
    opened = ndimage.binary_opening(sub, structure=np.ones((8, 8), bool))
    lab, n = ndimage.label(opened)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) >= 20:
            pts.append((float(cc.mean()) + 0.5, float(rr.mean()) + 0.5, int(len(rr))))

    # Also locate the 0 nM marker centered on the left spine y_col around y ~ 9.5..11.5
    c0 = int(round(fr["y_col"]))
    g = C.gray(image)
    r_lo0, r_hi0 = int(ya.to_pixel(11.6)), int(ya.to_pixel(9.2))
    # Left of y_col, marker protrudes outward (width >= 4 px to the left of left_band[0])
    lb0 = fr["left_band"][0]
    out_prof = (g[r_lo0:r_hi0, max(0, lb0 - 8) : lb0] < 130).sum(axis=1)
    hits0 = C.runs_1d(out_prof >= 3)
    if hits0:
        s0, e0 = max(hits0, key=lambda h: h[1] - h[0])
        r0_mid = r_lo0 + 0.5 * (s0 + e0) + 0.5
        pts.append((fr["y_col"], r0_mid, 100))

    pts.sort(key=lambda p: p[0])
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for c_m, r_m, area in pts:
        x_val = float(xa.to_value(c_m))
        if abs(x_val) < 0.3:
            x_val = 0.0
        y, lo, hi = C.val_interval(ya, r_m, 1.5 + cal_err)
        # Measure SD error bar (for x=0, look just right of spine if needed, or symmetric on c_m)
        if x_val == 0.0:
            # Cap protrudes left & right of left spine! Find top & bottom horizontal T-caps protruding left of left_band[0]
            lb0 = fr["left_band"][0]
            rb1 = fr["left_band"][1]
            # Top cap above marker (within 50 px), bottom cap below marker (within 50 px)
            up_caps = [r for r in range(int(r_m) - 55, int(r_m) - 8) if (g[r, max(0, lb0 - 7) : rb1 + 8] < 135).sum() >= 10]
            dn_caps = [r for r in range(int(r_m) + 9, min(int(rB), int(r_m) + 55)) if (g[r, max(0, lb0 - 7) : rb1 + 8] < 135).sum() >= 10]
            r_top_cap = float(np.median(up_caps)) + 0.5 if up_caps else (r_m - 30.0)
            r_bot_cap = float(np.median(dn_caps)) + 0.5 if dn_caps else (2.0 * r_m - r_top_cap)
            err_px = 1.5 + cal_err
            eb = {
                "bar_lo": float(ya.to_value(r_bot_cap)),
                "bar_hi": float(ya.to_value(r_top_cap)),
                "bar_lo_outer": float(ya.to_value(r_bot_cap + err_px)),
                "bar_hi_outer": float(ya.to_value(r_top_cap - err_px)),
                "px_lo": r_top_cap,
                "px_hi": r_bot_cap,
                "note": f"0 nM spine caps [{r_top_cap:.1f}, {r_bot_cap:.1f}]px",
            }
        else:
            eb = C.measure_error_bar_vertical(
                image, c_m, r_m, int(rT) + 4, int(rB) - 4, ya,
                dark_thresh=125.0, stem_half_w=3, gap_tol=3, mirror_lower_if_hidden=False,
            )
        out.append({
            "x": round(x_val, 3), "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_m, r_m], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"marker {area}px; {eb['note']}",
        })
    return out
