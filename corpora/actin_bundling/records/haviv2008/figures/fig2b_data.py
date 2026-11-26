"""Calibrated reader for Fig. 2b (p. 4 / 450): measured distribution Log(P(N_b)) = ln P(N_b) at phi = 3/1 (19 blue circles)."""
from __future__ import annotations

import math
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import NATIVE_QUANT_PX, calibrate_fig2b, value_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2b(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    blue = (B - R > 30.0) & (B - G > 30.0)
    mask = np.zeros_like(blue)
    mask[int(f["tr1"]) + 4 : int(f["br0"]) - 4, int(f["lc1"]) + 4 : int(f["rc1"]) + 45] = True
    blue &= mask

    lbl, _ = ndimage.label(blue)
    err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0
    readings: list[dict] = []
    for i, sl in enumerate(ndimage.find_objects(lbl), 1):
        sub = lbl[sl] == i
        area = int(sub.sum())
        if area < 400:
            continue
        cy, cx = ndimage.center_of_mass(sub)
        cx_f = float(cx + sl[1].start + 0.5)
        cy_f = float(cy + sl[0].start + 0.5)
        x_val = float(x_axis.to_value(cx_f))
        # Fig. 2b plots natural log ln P(N_b) linearly on y-axis; return linear P = exp(ln P)
        # so BoxGraph datum with scale: ln recovers exact ln P(N_b) and interval.
        ln_y, ln_lo, ln_hi = value_interval(y_axis, cy_f, err_px)
        readings.append({
            "x": x_val,
            "y": math.exp(ln_y),
            "lo": math.exp(ln_lo),
            "hi": math.exp(ln_hi),
            "px": [cx_f, cy_f],
            "px_lo": float(y_axis.to_pixel(ln_hi)),
            "px_hi": float(y_axis.to_pixel(ln_lo)),
            "note": f"blue circle area={area} px, lnP={ln_y:.4f} [{ln_lo:.4f}, {ln_hi:.4f}]; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    readings.sort(key=lambda r: r["x"])
    return readings
