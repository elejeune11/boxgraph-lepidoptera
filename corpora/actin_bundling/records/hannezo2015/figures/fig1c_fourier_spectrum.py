"""Calibrated reader for Main Fig. 1C bottom-right (p. 2): Fourier spectrum of larval actin rings."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_inward_ticks, match_ticks_regular, verify_axis, sample_solid_curve_y


def read(image, at):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=140.0)
    # Exclude top micrograph edge if present: take the two main plot-box horizontal lines
    h_box = [h for h in hl if h[0] > 50.0]
    h_top, h_bot = h_box[0], h_box[-1]
    # Take the two plot-box vertical lines (exclude panel outer border at right edge if present)
    v_box = [v for v in vl if v[0] < image.shape[1] - 12]
    v_left, v_right = v_box[0], v_box[-1]

    _, xt_bot, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 14), margin=10, max_width=10)
    # X ticks: left border = 0, interior ticks = 0.5, 1.0, 1.5, 2.0, 2.5, right border = 3.0
    xt_int = match_ticks_regular([t for t in xt_bot if t < v_right[0] - 20], 5, tol=3.0)
    x_marks = [(v_left[0], 0.0)] + list(zip(xt_int, [0.5, 1.0, 1.5, 2.0, 2.5])) + [(v_right[0], 3.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=1.5)

    # Y ticks: 8 interior ticks from 0.14 down to 0.00 (step 0.02), with 0.10 at 3rd tick from top and 0.0 at bottommost interior tick
    yt_int = match_ticks_regular([t for t in yt_left if t > h_top[0] + 20], 8, tol=2.5)
    y_vals = [0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02, 0.00]
    y_axis = verify_axis("linear", list(zip(yt_int, y_vals)), max_residual_px=1.5)

    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180) & (rgb[:, :, 1] < 90) & (rgb[:, :, 2] < 90)
    # Mask out legend line in top right
    red[: int(h_top[2] + 150), int(x_axis.to_pixel(1.8)) :] = False

    r_min, r_max = int(h_top[2]) + 10, int(h_bot[1]) - 5
    return [sample_solid_curve_y(red, x_axis, y_axis, x, r_min, r_max, search_dx=2) for x in at]
