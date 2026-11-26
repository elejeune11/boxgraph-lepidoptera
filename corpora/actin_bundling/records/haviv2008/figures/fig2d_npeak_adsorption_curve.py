"""Calibrated reader for Fig. 2d (p. 4 / 450): dashed black Langmuir-adsorption N_peak scaling curve."""
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

    dark = gray < 165.0
    for mc in marker_cols:
        dark[:, max(0, int(round(mc)) - 26) : min(image.shape[1], int(round(mc)) + 27)] = False

    r_lo = int(f["tr1"]) + 10
    r_hi = int(f["br0"]) - 30

    # Purely topological branch selection among Fig. 2d's three black curves:
    # - For 1.35 <= phi <= 2.42: dashed curve is the HIGHEST (smallest row index) curve in the panel.
    # - For 0.60 <= phi <= 0.90: dashed curve is the MIDDLE curve (strictly between dotted above and solid below).
    x_cross_dot = 1.25

    def _sample_col(c: int, x: float) -> tuple[float, float] | None:
        runs = col_runs(dark[r_lo:r_hi, c], r_lo)
        runs = [r for r in runs if (r[3] - r[2] + 1) <= 60]
        if not runs:
            return None
        runs.sort(key=lambda t: t[0])  # ascending row = top to bottom
        if x >= x_cross_dot:
            chosen = runs[0]  # highest curve
            # Verify it sits strictly above the lower curves (at least 2 runs in column or upper half)
            if len(runs) == 1 and chosen[0] > (r_lo + r_hi) / 2.0:
                return None
        else:
            if len(runs) >= 3:
                chosen = runs[1]  # middle of 3 curves
            elif len(runs) == 2:
                # If dotted curve is in a tiny dot gap at exact column c, verify a dot exists above runs[0] nearby
                above_dot = dark[r_lo : max(r_lo, int(runs[0][2]) - 40), max(int(f["lc1"]) + 10, c - 18) : min(int(f["rc0"]) - 10, c + 19)]
                if above_dot.any():
                    chosen = runs[0]
                else:
                    return None
            else:
                return None
        r0, r1 = chosen[2], chosen[3]
        w = 255.0 - gray[r0 : r1 + 1, c]
        rc = float(np.sum(np.arange(r0, r1 + 1) * w) / np.sum(w)) + 0.5
        return rc, chosen[1]

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        hit = _sample_col(c, float(x))
        used_c = c
        if hit is None:
            for dc in (2, -2, 4, -4, 6, -6):
                hit = _sample_col(c + dc, float(x))
                if hit is not None:
                    used_c = c + dc
                    break
        if hit is None:
            readings.append({"x": float(x), "y": None, "note": "in dash gap of dashed black curve"})
            continue
        rc, half_h = hit
        hl = _sample_col(used_c - 10, float(x))
        hr = _sample_col(used_c + 10, float(x))
        slope = abs(hr[0] - hl[0]) / 20.0 if (hl is not None and hr is not None) else 0.0
        err_px = half_h + NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + slope * (NATIVE_QUANT_PX + x_axis.residual_px + abs(used_c - c))
        y_val, y_lo, y_hi = value_interval(y_axis, rc, err_px)
        readings.append({
            "x": float(x),
            "y": y_val,
            "lo": max(0.0, y_lo),
            "hi": y_hi,
            "px": [float(used_c + 0.5), rc],
            "px_lo": float(y_axis.to_pixel(y_hi)),
            "px_hi": float(y_axis.to_pixel(max(0.0, y_lo))),
            "note": f"black dashed curve N_peak={y_val:.2f}, half-h={half_h:.1f} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
