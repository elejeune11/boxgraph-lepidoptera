"""Calibrated reader for Fig. 3(b) (p. 3): predicted overtwist vs (beta*U*rho)^2 for xi_t = 10 (olive)."""
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
    # Olive/yellow-green curve (xi_t = 10), median RGB ~= [168, 154, 87]
    olive = (r > b + 20.0) & (g > b + 20.0) & (r > 100.0) & (r < 225.0)
    return sample_curve_from_mask(
        olive, x_axis, y_axis, at, r_lo=1400, r_hi=r_ax - 2,
        search_dx=20, note_prefix="xi_t=10 olive; ",
    )
