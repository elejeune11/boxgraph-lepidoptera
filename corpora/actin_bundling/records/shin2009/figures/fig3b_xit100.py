"""Calibrated reader for Fig. 3(b) (p. 3): predicted overtwist vs (beta*U*rho)^2 for xi_t = 100 (purple)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import calibrate_fig3b, sample_curve_from_mask


def read(image, at):
    x_axis, y_axis, (c_ax, r_ax) = calibrate_fig3b(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    # Purple/magenta curve (xi_t = 100), median RGB ~= [157, 89, 131]
    purple = (r > g + 15.0) & (b > g + 10.0) & (np.abs(r - b) < 55.0) & (r > 100.0)
    return sample_curve_from_mask(
        purple, x_axis, y_axis, at, r_lo=120, r_hi=r_ax - 2,
        search_dx=20, note_prefix="xi_t=100 purple; ",
    )
