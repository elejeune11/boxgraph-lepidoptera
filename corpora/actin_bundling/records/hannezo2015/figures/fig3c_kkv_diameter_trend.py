"""Calibrated reader for Main Fig. 3C top-right (p. 4): smooth blue trend curve of tube diameter vs actin concentration."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_solid_curve_y
from fig3c_kkv_diameter_vs_actin import calibrate_fig3c_diam


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_fig3c_diam(image)
    rgb = image.astype(float)
    # Purple-blue guide curve
    blue = (rgb[:, :, 2] > 90) & (rgb[:, :, 2] - rgb[:, :, 1] > 30) & (rgb[:, :, 2] - rgb[:, :, 0] > 15)
    r_min, r_max = int(rt) + 10, int(rb) - 8
    return [sample_solid_curve_y(blue, x_axis, y_axis, x, r_min, r_max, search_dx=2) for x in at]
