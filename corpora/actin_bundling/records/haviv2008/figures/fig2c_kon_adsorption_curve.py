"""Calibrated reader for Fig. 2c (p. 4 / 450): red dashed Langmuir adsorption on-rate curve (2*k_on/D propto n_af)."""
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

    hline_r, hline_r0, hline_r1 = find_hline_row_red(image, f)
    marker_cols = detect_fig2c_marker_columns(image, f)

    red = (R - G > 20.0) & (R - B > 20.0)
    # Remove horizontal solid red line and marker/error-bar columns
    red[max(0, hline_r0 - 12) : min(image.shape[0], hline_r1 + 13), :] = False
    for mc in marker_cols:
        red[:, max(0, int(round(mc)) - 26) : min(image.shape[1], int(round(mc)) + 27)] = False

    r_lo = int(f["tr1"]) + 10
    r_hi = int(f["br0"]) - 30

    # Topological branch selection purely from the two non-horizontal red curves in Fig. 2c:
    # - Right of the crossing near phi ~ 1.15: the rising/saturating red dashed curve is the
    #   HIGHEST (smallest row index) red curve above the horizontal line.
    # - Left of the crossing (phi < 1.1): the rising red dashed curve is the LOWER (larger row index)
    #   red curve above the horizontal line (below the higher-peaking red dotted curve).
    x_cross = 1.15

    def _sample_branch(c: int, x: float) -> tuple[float, float] | None:
        # Restrict above horizontal line (all clean dashes for phi >= 0.4 lie above hline_r)
        runs = col_runs(red[r_lo : hline_r0 - 12, c], r_lo)
        # Exclude wide vertical artifacts
        runs = [r for r in runs if (r[3] - r[2] + 1) <= 45]
        if not runs:
            return None
        if x >= x_cross:
            chosen = min(runs, key=lambda t: t[0])  # topmost (highest y)
        else:
            # Left of crossing: if 2 branches exist above hline, pick the lower one (larger row)
            # If only 1 exists (because the dotted curve has a dot gap in this exact column),
            # verify it's part of a horizontal/sloped dash (width >= 16 px in local window)
            if len(runs) >= 2:
                chosen = max(runs, key=lambda t: t[0])
            else:
                # Check if another red run exists slightly above within +/-16 cols (the dotted curve)
                rc_cand = runs[0][0]
                above_dots = red[r_lo : max(r_lo, int(rc_cand) - 40), max(0, c - 18) : min(image.shape[1], c + 19)]
                if above_dots.any():
                    chosen = runs[0]
                else:
                    return None
        r0, r1 = chosen[2], chosen[3]
        w = (R[r0 : r1 + 1, c] - G[r0 : r1 + 1, c]) + (R[r0 : r1 + 1, c] - B[r0 : r1 + 1, c])
        rc = float(np.sum(np.arange(r0, r1 + 1) * w) / np.sum(w)) + 0.5
        return rc, chosen[1]

    readings: list[dict] = []
    for x in at:
        c = int(round(x_axis.to_pixel(x)))
        hit = _sample_branch(c, float(x))
        used_c = c
        if hit is None:
            # Check narrow +/-6 px tolerance in case x sits right at a dash edge
            for dc in (2, -2, 4, -4, 6, -6):
                hit = _sample_branch(c + dc, float(x))
                if hit is not None:
                    used_c = c + dc
                    break
        if hit is None:
            readings.append({"x": float(x), "y": None, "note": "in dash gap of red dashed curve"})
            continue
        rc, half_h = hit
        h_l = _sample_branch(used_c - 12, float(x))
        h_r = _sample_branch(used_c + 12, float(x))
        slope = abs(h_r[0] - h_l[0]) / 24.0 if (h_l is not None and h_r is not None) else 0.0
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
            "note": f"red dashed curve (2*kon/D={y_val:.4f}), half-h={half_h:.1f} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
