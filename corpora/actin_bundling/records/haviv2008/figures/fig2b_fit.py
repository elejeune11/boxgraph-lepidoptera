"""Calibrated reader for Fig. 2b (p. 4 / 450): Eq. (1) solid black fit curve on semi-log Log(P(N_b)) = ln P(N_b) plot."""
from __future__ import annotations

import math
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import NATIVE_QUANT_PX, calibrate_fig2b, col_runs, value_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2b(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0

    blue = (B - R > 22.0) & (B - G > 22.0)
    blue_dil = ndimage.binary_dilation(blue, iterations=4)

    r_lo = int(f["tr1"]) + 10
    r_hi = int(f["br0"]) - 30  # exclude bottom inward ticks
    black = (gray < 150.0) & (B - R < 20.0) & (B - G < 20.0) & (~blue_dil)

    def _sample_col(c: int) -> tuple[float, float] | None:
        runs = col_runs(black[r_lo:r_hi, c], r_lo)
        if not runs:
            return None
        best = max(runs, key=lambda t: t[1])
        r0, r1 = best[2], best[3]
        w = 255.0 - gray[r0 : r1 + 1, c]
        rc = float(np.sum(np.arange(r0, r1 + 1) * w) / np.sum(w)) + 0.5
        return rc, best[1]

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        hit = _sample_col(c)
        if hit is None:
            readings.append({"x": float(x), "y": None, "note": "black curve occluded at this N_b"})
            continue
        rc, half_h = hit
        h_l = _sample_col(max(int(f["lc1"]) + 12, c - 16))
        h_r = _sample_col(min(int(f["rc0"]) - 12, c + 16))
        slope = abs(h_r[0] - h_l[0]) / 32.0 if (h_l is not None and h_r is not None) else 0.0
        err_px = half_h + NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + slope * (NATIVE_QUANT_PX + x_axis.residual_px)
        ln_y, ln_lo, ln_hi = value_interval(y_axis, rc, err_px)
        readings.append({
            "x": float(x),
            "y": math.exp(ln_y),
            "lo": math.exp(ln_lo),
            "hi": math.exp(ln_hi),
            "px": [float(c + 0.5), rc],
            "px_lo": float(y_axis.to_pixel(ln_hi)),
            "px_hi": float(y_axis.to_pixel(ln_lo)),
            "note": f"black curve lnP={ln_y:.4f} [{ln_lo:.4f}, {ln_hi:.4f}], half-h={half_h:.1f} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
