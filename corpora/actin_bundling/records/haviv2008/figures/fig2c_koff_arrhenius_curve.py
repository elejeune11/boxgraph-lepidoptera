"""Calibrated reader for Fig. 2c (p. 4 / 450): black solid Arrhenius off-rate fit curve (10*k_off/D divided by 10)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import (
    NATIVE_QUANT_PX,
    calibrate_fig2c,
    col_runs,
    detect_fig2c_marker_columns,
    find_hline_row_red,
    value_interval,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2c(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0

    hline_r, hline_r0, hline_r1 = find_hline_row_red(image, f)
    marker_cols = detect_fig2c_marker_columns(image, f)

    r_lo = int(f["tr1"]) + 10
    r_hi = int(f["br0"]) - 30  # stay clear of bottom inward ticks

    # Mask out marker/error-bar columns and horizontal red line
    clean_cols = np.ones(image.shape[1], dtype=bool)
    for mc in marker_cols:
        clean_cols[max(0, int(round(mc)) - 26) : min(image.shape[1], int(round(mc)) + 27)] = False

    # Where the black curve lies ABOVE the red horizontal line (left part, phi < 2.5),
    # it is strictly neutral black and separated from all red curves.
    # Where it lies BELOW the red horizontal line (right part, phi > 2.5), it is the ONLY dark
    # curve below hline_r (with the red dotted curve superimposed along its upper edge).
    neutral_black = (gray < 150.0) & (R - G < 25.0) & (R - B < 25.0)
    dark_below = gray < 155.0
    dark_below[: hline_r1 + 14, :] = False

    def _sample_col(c: int) -> tuple[float, float, bool] | None:
        if not clean_cols[c]:
            return None
        # Check neutral black above horizontal line first
        above_runs = col_runs(neutral_black[r_lo : hline_r0 - 10, c], r_lo)
        below_runs = col_runs(dark_below[hline_r1 + 14 : r_hi, c], hline_r1 + 14)
        if above_runs:
            best = max(above_runs, key=lambda t: t[1])
            r0, r1 = best[2], best[3]
            w = 255.0 - gray[r0 : r1 + 1, c]
            rc = float(np.sum(np.arange(r0, r1 + 1) * w) / np.sum(w)) + 0.5
            return rc, best[1], False
        if below_runs:
            best = max(below_runs, key=lambda t: t[1])
            r0, r1 = best[2], best[3]
            w = 255.0 - gray[r0 : r1 + 1, c]
            rc = float(np.sum(np.arange(r0, r1 + 1) * w) / np.sum(w)) + 0.5
            # Check if red dotted curve overlaps here
            has_red = bool(((R[r0 : r1 + 1, c] - G[r0 : r1 + 1, c]) > 22.0).any())
            return rc, best[1], has_red
        # Near the crossing with the red horizontal line (y ~ 0.45): check neutral black across full height
        all_nb = col_runs(neutral_black[r_lo:r_hi, c], r_lo)
        if all_nb:
            best = max(all_nb, key=lambda t: t[1])
            r0, r1 = best[2], best[3]
            w = 255.0 - gray[r0 : r1 + 1, c]
            rc = float(np.sum(np.arange(r0, r1 + 1) * w) / np.sum(w)) + 0.5
            return rc, best[1], False
        return None

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        hit = _sample_col(c)
        if hit is None:
            readings.append({"x": float(x), "y": None, "note": "black Arrhenius curve not found at this phi"})
            continue
        rc, half_h, overlap_red = hit
        h_l = _sample_col(max(int(f["lc1"]) + 10, c - 16))
        h_r = _sample_col(min(int(f["rc0"]) - 10, c + 16))
        slope = abs(h_r[0] - h_l[0]) / 32.0 if (h_l is not None and h_r is not None) else 0.0
        extra = 4.0 if overlap_red else 0.0
        err_px = half_h + NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + slope * (NATIVE_QUANT_PX + x_axis.residual_px) + extra
        y10, y10_lo, y10_hi = value_interval(y_axis, rc, err_px)
        readings.append({
            "x": float(x),
            "y": y10 / 10.0,
            "lo": max(0.0, y10_lo / 10.0),
            "hi": y10_hi / 10.0,
            "px": [float(c + 0.5), rc],
            "px_lo": float(y_axis.to_pixel(y10_hi)),
            "px_hi": float(y_axis.to_pixel(max(0.0, y10_lo))),
            "note": f"black solid curve (10*koff/D={y10:.4f}, div by 10){'; red dotted curve coalesced' if overlap_red else ''}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
