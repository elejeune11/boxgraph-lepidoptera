"""Calibrated reader for Main Fig. 2A bottom-right (p. 3): FRAP linear variance Model line."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_dashed_curve_y
from fig2a_frap_variance import calibrate_fig2a_diff


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_fig2a_diff(image)
    rgb = image.astype(float)
    green = (rgb[:, :, 1] > 120) & (rgb[:, :, 1] - rgb[:, :, 2] > 50) & (rgb[:, :, 0] < 170)
    green[: int(rt) + 280, int(cl) + int(0.35 * (cr - cl)) :] = False  # exclude legend

    r_min, r_max = int(rt) + 12, int(rb) - 4
    return [sample_dashed_curve_y(green, x_axis, y_axis, x, r_min, r_max, max_search_dx=35) for x in at]
