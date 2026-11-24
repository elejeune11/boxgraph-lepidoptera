"""Shared calibrated reader for Fig. 6B stacked abundance bar chart (1200 dpi)."""

from __future__ import annotations
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

K_UNBIND_VALUES = [
    0.1, 1.0, 1.125, 1.15625, 1.1875, 1.21875, 1.25,
    1.28125, 1.3125, 1.34375, 1.375, 1.5, 1.625, 1.75,
]


def read_fig6b(image: np.ndarray, category: str) -> list[dict]:
    """Read one stacked bar category ('ring', 'shell_sb', or 'shell') across all 14 k_unbind bars."""
    info = pc.find_plot_box(image)
    y_axis = pc.fit_axis(info["y_ticks"], [100.0, 75.0, 50.0, 25.0, 0.0])
    x_ticks = info["x_ticks"]
    if len(x_ticks) != len(K_UNBIND_VALUES):
        raise ValueError(f"Expected 14 x-ticks in Fig. 6B, got {len(x_ticks)}")
    c0, r0, c1, r1 = info["box"]
    px_per_pct = abs(y_axis.to_pixel(0.0) - y_axis.to_pixel(100.0)) / 100.0
    y_misfit = max(y_axis.residual_px, 0.5)
    # Border lines inside stacked bars are ~3 px thick; half-uncertainty covers border + calibration
    err_pct = (2.5 + y_misfit) / px_per_pct

    readings = []
    for k_val, xt in zip(K_UNBIND_VALUES, x_ticks):
        cx = int(round(xt))
        # Use central 24 columns of each bar (bar width ~82 px), staying well clear of left/right black borders
        strip = image[r0 + 2:r1 - 1, cx - 12:cx + 12].astype(np.float64)
        col_rgb = strip.mean(axis=1)
        is_green = (col_rgb[:, 1] - col_rgb[:, 0] > 25) & (col_rgb[:, 1] > 110)
        is_yellow = (col_rgb[:, 0] > 165) & (col_rgb[:, 1] > 135) & (col_rgb[:, 2] < 85)
        is_red = (col_rgb[:, 0] > 115) & (col_rgb[:, 1] < 75) & (col_rgb[:, 2] < 75)

        r_red = np.flatnonzero(is_red) + (r0 + 2)
        r_yel = np.flatnonzero(is_yellow) + (r0 + 2)
        r_grn = np.flatnonzero(is_green) + (r0 + 2)

        # Compute exact transition rows (stack order bottom->top: red Shell, yellow Shell-SB, green Ring)
        row_0 = y_axis.to_pixel(0.0)
        row_100 = y_axis.to_pixel(100.0)

        # Shell (red) occupies [row_red_top, row_0]
        if r_red.size >= 8:
            if r_yel.size >= 8:
                row_red_top = 0.5 * (float(r_red.min()) + float(r_yel.max()) + 1.0)
            elif r_grn.size >= 8:
                row_red_top = 0.5 * (float(r_red.min()) + float(r_grn.max()) + 1.0)
            else:
                row_red_top = row_100
            y_shell = min(100.0, max(0.0, y_axis.to_value(row_red_top)))
        else:
            row_red_top = row_0
            y_shell = 0.0

        # Ring (green) occupies [row_100, row_grn_bot]
        if r_grn.size >= 8:
            if r_yel.size >= 8:
                row_grn_bot = 0.5 * (float(r_grn.max()) + 1.0 + float(r_yel.min()))
            elif r_red.size >= 8:
                row_grn_bot = 0.5 * (float(r_grn.max()) + 1.0 + float(r_red.min()))
            else:
                row_grn_bot = row_0
            y_ring = min(100.0, max(0.0, 100.0 - y_axis.to_value(row_grn_bot)))
        else:
            row_grn_bot = row_100
            y_ring = 0.0

        # Shell w. strong bundling (yellow) occupies the middle/remaining band
        if r_yel.size >= 8:
            y_shell_sb = min(100.0, max(0.0, 100.0 - y_shell - y_ring))
        else:
            y_shell_sb = 0.0

        if category == "shell":
            val = y_shell
            px_row = row_red_top if val > 0 else row_0
            note = f"Fig. 6B Shell (red) bar height (ry={y_axis.residual_px:.2f}px)"
        elif category == "shell_sb":
            val = y_shell_sb
            # Place marker at top edge of yellow band (which is row_grn_bot if ring>0, else row_100 if yel>0, else row_0)
            px_row = row_grn_bot if val > 0 else row_0
            note = f"Fig. 6B Shell-SB (yellow) stacked height; top at {y_shell + val:.1f}% (ry={y_axis.residual_px:.2f}px)"
        elif category == "ring":
            val = y_ring
            # Place marker at bottom edge of green Ring band (transition from yellow/red to green, or top spine if 0%)
            px_row = row_grn_bot if val > 0 else row_100
            note = f"Fig. 6B Ring (green) stacked height; starts at {100.0 - val:.1f}% (ry={y_axis.residual_px:.2f}px)"
        else:
            raise ValueError(f"Unknown category {category!r}")

        lo = max(0.0, val - err_pct)
        hi = min(100.0, val + err_pct)
        h_px = err_pct * px_per_pct
        readings.append({
            "x": float(k_val),
            "y": float(val),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(xt), float(px_row)],
            "px_lo": float(px_row - h_px),
            "px_hi": float(px_row + h_px),
            "note": note,
        })
    return readings
