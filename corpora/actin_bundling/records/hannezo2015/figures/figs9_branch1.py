"""Calibrated reader for SI Fig. S9 (p. 16): primary stable branch (i = 1, solid blue) chi/phi vs vertical axis."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_horizontal_x_of_y
from figs9_bifurcation_points import calibrate_figs9


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs9(image)
    rgb = image.astype(float)
    blue = (rgb[:, :, 2] > 180.0) & (rgb[:, :, 0] < 65.0) & (rgb[:, :, 1] < 65.0)

    # Mask out the 4 density-profile inset boxes so inset blue curves never interfere:
    # Branch 1 lies strictly in chi/phi in [1.05, 1.55]; inset 1 lies left of chi/phi = 1.12 and y in [3.2, 4.7]
    c_ins1 = int(x_axis.to_pixel(1.14))
    r_ins1_top = int(y_axis.to_pixel(4.75))
    r_ins1_bot = int(y_axis.to_pixel(3.15))
    blue[r_ins1_top:r_ins1_bot, :c_ins1] = False

    c_min = int(x_axis.to_pixel(1.05))
    c_max = int(x_axis.to_pixel(1.48))

    # Pick leftmost blue run in [1.05, 1.48] (branch i = 1 is solid and well-separated from branch i = 2)
    def pick_leftmost(runs, r):
        return min(runs, key=lambda t: t[0])

    return [sample_horizontal_x_of_y(blue, x_axis, y_axis, y_val, c_min, c_max, pick_fn=pick_leftmost) for y_val in at]
