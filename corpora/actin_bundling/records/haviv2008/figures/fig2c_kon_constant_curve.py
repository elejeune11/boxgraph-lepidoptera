"""Calibrated reader for Fig. 2c (p. 4 / 450): red solid horizontal constant on-rate line (2*k_on/D propto n_a)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import NATIVE_QUANT_PX, calibrate_fig2c, find_hline_row_red, value_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2c(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    red = (R - G > 20.0) & (R - B > 20.0)

    hline_r, hline_r0, hline_r1 = find_hline_row_red(image, f)
    half_h_global = (hline_r1 - hline_r0 + 1) / 2.0

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        r_s = max(int(f["tr1"]) + 5, hline_r0 - 8)
        r_e = min(int(f["br0"]) - 5, hline_r1 + 9)
        col_m = red[r_s:r_e, c]
        idx = np.where(col_m)[0]
        if len(idx) == 0:
            rc = hline_r
            half_h = half_h_global
        else:
            rows = idx + r_s
            w = (R[rows, c] - G[rows, c]) + (R[rows, c] - B[rows, c])
            rc = float(np.sum(rows * w) / np.sum(w)) + 0.5
            half_h = (int(rows[-1]) - int(rows[0]) + 1) / 2.0
        err_px = half_h + NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0)
        y_val, y_lo, y_hi = value_interval(y_axis, rc, err_px)
        readings.append({
            "x": float(x),
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "px": [float(c + 0.5), rc],
            "px_lo": float(y_axis.to_pixel(y_hi)),
            "px_hi": float(y_axis.to_pixel(y_lo)),
            "note": f"red solid horizontal line (2*kon/D={y_val:.4f}), half-h={half_h:.1f} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
