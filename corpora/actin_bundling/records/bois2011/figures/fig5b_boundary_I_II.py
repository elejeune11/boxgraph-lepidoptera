"""Calibrated reader for Fig. 5(b): ASDM Region I/II stability boundary (Pe = 7, rho_a = 1)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy import ndimage

from _calib import (
    find_frame,
    find_inward_ticks,
    sample_horizontal_x_of_y,
    verify_axis,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read critical Da/Ds horizontally at each vertical-axis rho_s > 1 in `at`."""
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45, thresh=110.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(dark, h_top, h_bot, v_left, v_right)
    mid_x = 0.5 * (v_left[0] + v_right[0])
    x1 = min(xt_bot, key=lambda p: abs(p - mid_x))
    x_axis = verify_axis("linear", [(v_left[0], 0.0), (x1, 1.0), (v_right[0], 2.0)], max_residual_px=1.0)

    mid_y = 0.5 * (h_top[0] + h_bot[0])
    y1 = min(yt_left, key=lambda p: abs(p - mid_y))
    y_axis = verify_axis("linear", [(h_top[0], 2.0), (y1, 1.0), (h_bot[0], 0.0)], max_residual_px=1.0)

    # Isolate the continuous Region I/II curve in rho_s in (1, 2], Da/Ds in [0.7, 1.5]
    r_lo = h_top[2] + 6
    r_hi = int(round(y_axis.to_pixel(1.0))) - 12
    c_lo = int(round(x_axis.to_pixel(0.70)))
    c_hi = int(round(x_axis.to_pixel(1.48)))
    sub = np.zeros_like(dark)
    sub[r_lo : r_hi + 1, c_lo : c_hi + 1] = dark[r_lo : r_hi + 1, c_lo : c_hi + 1]

    lab, n = ndimage.label(sub)
    counts = np.bincount(lab.ravel())
    counts[0] = 0
    best_id = int(np.argmax(counts))
    curve_mask = lab == best_id

    out = []
    for rs in at:
        res = sample_horizontal_x_of_y(
            curve_mask,
            x_axis,
            y_axis,
            float(rs),
            c_lo,
            c_hi,
            note_prefix="Fig. 5(b) I/II boundary, ",
        )
        out.append(res)
    return out
