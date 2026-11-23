"""Calibrated reader for SI Fig. S12 (p. 20): primary supercritical branch (i_c = 6, alpha = 15, phi = 3), chi vs max_rho."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_horizontal_x_of_y
from figs12_bifurcation_points import calibrate_figs12


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs12(image)
    rgb = image.astype(float)
    ink = rgb.min(axis=2) < 180.0
    c_min, c_max = int(cl) + 25, int(cr) - 25

    # The primary supercritical i_c = 6 branch is the leftmost curve across max_rho in [1.3, 13.5]
    def pick_leftmost(runs, r):
        return runs[0] if runs else None

    return [sample_horizontal_x_of_y(ink, x_axis, y_axis, y_val, c_min, c_max, pick_fn=pick_leftmost) for y_val in at]
