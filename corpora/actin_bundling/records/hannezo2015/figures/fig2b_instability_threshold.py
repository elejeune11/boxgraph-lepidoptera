"""Calibrated reader for Main Fig. 2B (p. 3): phase diagram instability threshold chi_c vs turnover phi."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_inward_ticks, match_ticks_regular, verify_axis, sample_solid_curve_y


def read(image, at):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=140.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 15), margin=12, max_width=10)

    # X ticks: left border = 0, interior ticks = 2, 4, 6, 8, right border = 10
    xt_int = match_ticks_regular(xt_top, 4, tol=3.0)
    x_axis = verify_axis("linear", [(v_left[0], 0.0)] + list(zip(xt_int, [2.0, 4.0, 6.0, 8.0])) + [(v_right[0], 10.0)], max_residual_px=1.5)

    # Y ticks: top border = 20, interior ticks = 15, 10, 5, bottom border = 0
    yt_int = match_ticks_regular(yt_left, 3, tol=3.0)
    y_axis = verify_axis("linear", [(h_top[0], 20.0)] + list(zip(yt_int, [15.0, 10.0, 5.0])) + [(h_bot[0], 0.0)], max_residual_px=1.5)

    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180) & (rgb[:, :, 1] < 60) & (rgb[:, :, 2] < 60)
    # Exclude legend line in top right
    red[: int(h_top[2] + 150), int(x_axis.to_pixel(7.0)) :] = False

    r_min, r_max = int(h_top[2]) + 10, int(h_bot[1]) - 5
    return [sample_solid_curve_y(red, x_axis, y_axis, x, r_min, r_max) for x in at]
