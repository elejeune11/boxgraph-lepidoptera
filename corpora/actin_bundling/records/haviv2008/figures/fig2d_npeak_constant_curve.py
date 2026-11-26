"""Calibrated reader for Fig. 2d (p. 4 / 450): solid black constant-on-rate N_peak scaling curve."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import (
    NATIVE_QUANT_PX,
    calibrate_fig2d,
    col_runs,
    detect_fig2d_marker_columns,
    value_interval,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2d(image)
    gray = image.astype(float).mean(axis=2)
    marker_cols = detect_fig2d_marker_columns(image, f)

    dark = gray < 160.0
    for mc in marker_cols:
        dark[:, max(0, int(round(mc)) - 26) : min(image.shape[1], int(round(mc)) + 27)] = False

    r_lo = int(f["tr1"]) + 10
    r_hi = int(f["br0"]) - 30  # stay clear of bottom inward ticks

    # Purely topological branch selection among Fig. 2d's three black curves:
    # - For 0.5 < phi < 2.2: solid curve is the LOWEST (largest row index) curve in the panel.
    # - For 2.55 < phi < 3.55: solid curve is the HIGHEST (smallest row index) curve (above dotted,
    #   after the steeper dashed curve has already exited the top frame around phi ~ 2.45).
    x_cross_dot = 2.40

    def _sample_col(c: int, x: float) -> tuple[float, float] | None:
        runs = col_runs(dark[r_lo:r_hi, c], r_lo)
        runs = [r for r in runs if (r[3] - r[2] + 1) <= 160]
        if not runs:
            return None
        chosen = max(runs, key=lambda t: t[0]) if x < x_cross_dot else min(runs, key=lambda t: t[0])
        r0, r1 = chosen[2], chosen[3]
        w = 255.0 - gray[r0 : r1 + 1, c]
        rc = float(np.sum(np.arange(r0, r1 + 1) * w) / np.sum(w)) + 0.5
        return rc, chosen[1]

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        hit = _sample_col(c, float(x))
        if hit is None:
            readings.append({"x": float(x), "y": None, "note": "solid curve not found at this phi"})
            continue
        rc, half_h = hit
        hl = _sample_col(c - 12, float(x))
        hr = _sample_col(c + 12, float(x))
        slope = abs(hr[0] - hl[0]) / 24.0 if (hl is not None and hr is not None) else 0.0
        err_px = half_h + NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + slope * (NATIVE_QUANT_PX + x_axis.residual_px)
        y_val, y_lo, y_hi = value_interval(y_axis, rc, err_px)
        readings.append({
            "x": float(x),
            "y": y_val,
            "lo": max(0.0, y_lo),
            "hi": y_hi,
            "px": [float(c + 0.5), rc],
            "px_lo": float(y_axis.to_pixel(y_hi)),
            "px_hi": float(y_axis.to_pixel(max(0.0, y_lo))),
            "note": f"black solid curve N_peak={y_val:.2f}, half-h={half_h:.1f} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
