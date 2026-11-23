"""Calibrated reader for SI Fig. S11 (p. 19): primary subcritical branch (i_c = 6, alpha = 35, phi = 0.1), chi vs max_rho."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import sample_horizontal_x_of_y
from figs11_bifurcation_points import calibrate_figs11


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs11(image)
    rgb = image.astype(float)
    # Any non-white curve ink (curves change AUTO color across secondary bifurcations)
    ink = rgb.min(axis=2) < 180.0
    c_min, c_max = int(cl) + 25, int(cr) - 25

    out = []
    for y_val in at:
        # Below the secondary bifurcation at max_rho ~ 3.8, the primary i_c = 6 subcritical leg is the
        # leftmost curve; above max_rho ~ 6.0, a secondary symmetry-breaking branch lies slightly to its left,
        # so the primary i_c = 6 continuation is the 2nd curve from the left (red branch).
        if y_val < 3.6:
            def pick(runs, r):
                return runs[0] if runs else None
            prefix = "lower subcritical i=6 leg; "
        else:
            def pick(runs, r):
                return runs[1] if len(runs) >= 2 else None
            prefix = "upper i=6 branch (2nd from left); "
        out.append(sample_horizontal_x_of_y(ink, x_axis, y_axis, y_val, c_min, c_max, pick_fn=pick, note_prefix=prefix))
    return out
