"""Calibrated reader for Fig. 3(b): single-species linear stability boundary Pe/(1+(beta*pi*l/L)^2) vs c0."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from _calib import (
    find_frame,
    find_inward_ticks,
    match_ticks_regular,
    runs_1d,
    sample_horizontal_x_of_y,
    verify_axis,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read rescaled critical Peclet number horizontally at each vertical-axis c0 in `at`."""
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45, thresh=100.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    # Major ticks extend further inward (~16 px) than minor log ticks (~8 px); check band 9..15 px inside border
    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(
        dark, h_top, h_bot, v_left, v_right, band=(9, 15), margin=20, max_width=14, frac=0.6
    )
    # x axis (log10): left border (1), 1 interior major tick (10), right border (100)
    if len(xt_bot) != 1:
        mid = 0.5 * (v_left[0] + v_right[0])
        x10 = min(xt_bot, key=lambda p: abs(p - mid))
    else:
        x10 = xt_bot[0]
    x_axis = verify_axis("log10", [(v_left[0], 1.0), (x10, 10.0), (v_right[0], 100.0)], max_residual_px=1.0)

    # y axis (log10): top border (100), 3 interior major ticks (10, 1, 0.1), bottom border (0.01)
    y_inner = match_ticks_regular(yt_left, 3, tol=2.0)
    y_axis = verify_axis(
        "log10",
        [(h_top[0], 100.0), (y_inner[0], 10.0), (y_inner[1], 1.0), (y_inner[2], 0.1), (h_bot[0], 0.01)],
        max_residual_px=1.0,
    )

    c0 = v_left[2] + 15
    c1 = v_right[1] - 15
    # Build clean boundary mask: thick dark curve separating white left (gray > 245) from gray-shaded right (gray < 235)
    boundary = np.zeros_like(dark)
    H, W = dark.shape
    for r in range(h_top[2] + 10, h_bot[1] - 10):
        # Find dark runs in this row
        rs = [(a + c0, b + c0) for a, b in runs_1d(dark[r, c0 : c1 + 1])]
        for a, b in rs:
            if (b - a + 1) < 10:
                continue
            # Check white to the left (10..35 px left) and gray shading to the right (8..30 px right)
            l_lo, l_hi = max(c0, a - 35), max(c0 + 1, a - 8)
            r_lo, r_hi = min(c1 - 1, b + 8), min(c1, b + 35)
            if l_hi > l_lo and r_hi > r_lo:
                left_med = float(np.median(gray[r, l_lo:l_hi]))
                right_med = float(np.median(gray[r, r_lo:r_hi]))
                if left_med > 245.0 and 180.0 < right_med < 235.0:
                    boundary[r, a : b + 1] = True
                    break

    out = []
    for c0_val in at:
        res = sample_horizontal_x_of_y(
            boundary,
            x_axis,
            y_axis,
            float(c0_val),
            c0,
            c1,
            note_prefix="Fig. 3(b) boundary, ",
        )
        out.append(res)
    return out
