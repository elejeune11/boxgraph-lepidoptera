"""Calibrated reader for SI Fig. S9 (p. 16): second branch (i = 2, dashed blue) chi/phi vs vertical axis."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_horizontal_dashed_x_of_y
from figs9_bifurcation_points import calibrate_figs9


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs9(image)
    rgb = image.astype(float)
    blue = (rgb[:, :, 2] > 180.0) & (rgb[:, :, 0] < 65.0) & (rgb[:, :, 1] < 65.0)

    # Mask out inset 2 (chi/phi in [2.05, 2.65], y in [5.1, 6.4]) and inset 3 ([2.75, 3.35], y in [3.9, 5.2])
    blue[int(y_axis.to_pixel(6.45)) : int(y_axis.to_pixel(5.05)), int(x_axis.to_pixel(2.05)) : int(x_axis.to_pixel(2.68))] = False
    blue[int(y_axis.to_pixel(5.25)) : int(y_axis.to_pixel(3.85)), int(x_axis.to_pixel(2.75)) : int(x_axis.to_pixel(3.38))] = False

    # Branch 2 runs cleanly in a diagonal corridor between branch 1 and branch 3
    out = []
    for y_val in at:
        # Linear corridor centered around branch 2: chi/phi ~ 1.39 + 0.12*(y - 1)
        x_ctr = 1.39 + 0.122 * (y_val - 1.0)
        c_min = int(x_axis.to_pixel(max(1.36, x_ctr - 0.14)))
        c_max = int(x_axis.to_pixel(min(2.15, x_ctr + 0.14)))
        out.append(sample_horizontal_dashed_x_of_y(blue, x_axis, y_axis, y_val, c_min, c_max, max_search_dy=55))
    return out
