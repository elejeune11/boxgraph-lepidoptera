"""Calibrated reader for Fig. 3(b) (p. 3): predicted overtwist vs (beta*U*rho)^2 for xi_t = 1000 (blue)."""
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
    # Blue curve (xi_t = 1000), median RGB ~= [77, 88, 165]
    blue = (b > r + 15.0) & (b > g + 15.0) & (b > 100.0)
    return sample_curve_from_mask(
        blue, x_axis, y_axis, at, r_lo=120, r_hi=r_ax - 2,
        search_dx=20, note_prefix="xi_t=1000 blue; ",
    )
