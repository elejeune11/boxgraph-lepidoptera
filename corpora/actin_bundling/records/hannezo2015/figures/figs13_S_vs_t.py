"""Calibrated reader for SI Fig. S13 main plot (p. 22): order parameter S(t) vs time (s)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_inward_ticks, match_ticks_regular, verify_axis, sample_solid_curve_y


def calibrate_figs13_main(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=130.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    _, xt_bot, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(6, 22), margin=15, max_width=14)

    # X ticks: left=0, 5 interior = 500, 1000, 1500, 2000, 2500, right=3000
    xt_int = match_ticks_regular(xt_bot, 5, tol=3.5)
    x_axis = verify_axis(
        "linear",
        [(v_left[0], 0.0)] + list(zip(xt_int, [500.0, 1000.0, 1500.0, 2000.0, 2500.0])) + [(v_right[0], 3000.0)],
        max_residual_px=3.0,
    )
    # Y ticks: top=14000, 6 interior = 12000, 10000, 8000, 6000, 4000, 2000, bottom=0
    yt_int = match_ticks_regular(yt_left, 6, tol=4.0)
    y_axis = verify_axis(
        "linear",
        [(h_top[0], 14000.0)] + list(zip(yt_int, [12000.0, 10000.0, 8000.0, 6000.0, 4000.0, 2000.0])) + [(h_bot[0], 0.0)],
        max_residual_px=4.0,
    )
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs13_main(image)
    rgb = image.astype(float)
    blue = (rgb[:, :, 2] > 180.0) & (rgb[:, :, 0] < 60.0) & (rgb[:, :, 1] < 60.0)
    # Mask out the inset box (time in [1050, 2700], S in [2800, 9600])
    blue[int(y_axis.to_pixel(9700.0)) : int(y_axis.to_pixel(2600.0)), int(x_axis.to_pixel(1000.0)) : int(x_axis.to_pixel(2720.0))] = False

    r_min, r_max = int(rt) + 4, int(rb) - 2
    return [sample_solid_curve_y(blue, x_axis, y_axis, x, r_min, r_max, search_dx=3) for x in at]
