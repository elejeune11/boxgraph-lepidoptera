"""Calibrated reader for SI Fig. S15 (p. 22): red square-root fit curve of steady-state ring count vs chi."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_solid_curve_y
from figs15_n_rings_vs_chi import calibrate_figs15


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs15(image)
    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180.0) & (rgb[:, :, 1] < 80.0) & (rgb[:, :, 2] < 80.0)
    r_min, r_max = int(rt) + 6, int(rb) - 20
    return [sample_solid_curve_y(red, x_axis, y_axis, x, r_min, r_max, search_dx=3) for x in at]
