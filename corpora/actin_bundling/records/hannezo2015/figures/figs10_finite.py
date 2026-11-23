"""Calibrated reader for SI Fig. S10 (p. 17): finite-system number of patterns i_c/2 vs phi at alpha = 25."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_inward_ticks, match_ticks_regular, verify_axis, sample_solid_curve_y


def calibrate_figs10(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=130.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(6, 20), margin=15, max_width=12)
    xt_maj = match_ticks_regular(xt_top, 5, tol=3.0)
    x_axis = verify_axis("linear", list(zip(xt_maj, [0.0, 0.5, 1.0, 1.5, 2.0])), max_residual_px=1.5)
    yt_maj = match_ticks_regular(yt_left, 5, tol=4.0)
    y_axis = verify_axis("linear", list(zip(yt_maj, [4.0, 3.0, 2.0, 1.0, 0.0])), max_residual_px=3.0)
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs10(image)
    rgb = image.astype(float)
    # Blue step curve ("Finite system"): B - R > 20 and B > 100
    blue = (rgb[:, :, 2] - rgb[:, :, 0] > 20.0) & (rgb[:, :, 2] > 100.0)
    # Exclude legend box (phi in [1.0, 1.9], y in [2.0, 3.3])
    c_leg0 = int(x_axis.to_pixel(1.02))
    c_leg1 = int(x_axis.to_pixel(1.90))
    r_leg0 = int(y_axis.to_pixel(3.30))
    r_leg1 = int(y_axis.to_pixel(2.00))
    blue[r_leg0:r_leg1, c_leg0:c_leg1] = False

    r_min, r_max = int(rt) + 12, int(rb) - 12
    return [sample_solid_curve_y(blue, x_axis, y_axis, x, r_min, r_max, search_dx=2) for x in at]
