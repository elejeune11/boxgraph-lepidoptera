"""Calibrated reader for SI Fig. S10 (p. 17): infinite-system number of patterns alpha*phi^(1/4)/(2*pi) vs phi at alpha = 25."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_solid_curve_y
from figs10_finite import calibrate_figs10


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs10(image)
    rgb = image.astype(float)
    # Purple/maroon smooth curve ("Infinite system"): R - G > 30 and R >= B and R > 100
    purple = (rgb[:, :, 0] - rgb[:, :, 1] > 30.0) & (rgb[:, :, 0] >= rgb[:, :, 2]) & (rgb[:, :, 0] > 100.0)
    # Exclude legend box
    c_leg0 = int(x_axis.to_pixel(1.02))
    c_leg1 = int(x_axis.to_pixel(1.90))
    r_leg0 = int(y_axis.to_pixel(3.30))
    r_leg1 = int(y_axis.to_pixel(2.00))
    purple[r_leg0:r_leg1, c_leg0:c_leg1] = False

    r_min, r_max = int(rt) + 12, int(rb) - 12
    return [sample_solid_curve_y(purple, x_axis, y_axis, x, r_min, r_max, search_dx=2) for x in at]
