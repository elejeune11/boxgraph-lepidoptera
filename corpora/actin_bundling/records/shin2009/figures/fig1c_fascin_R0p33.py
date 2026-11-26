"""Calibrated reader for Fig. 1(c) (p. 2): fascin:actin SAXS intensity I(q) at R = 0.33 (magenta curve)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import calibrate_fig1c, sample_curve_from_mask


def read(image, at):
    x_axis, y_axis, (cL, rT, cR, rB) = calibrate_fig1c(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    # Magenta/purple curve (R = 0.33)
    mag = (r - g > 35.0) & (b - g > 30.0) & (r > 110.0) & (b > 110.0)
    return sample_curve_from_mask(
        mag, x_axis, y_axis, at, r_lo=rT + 10, r_hi=rB - 10,
        search_dx=28, note_prefix="fig1c fascin R=0.33 magenta; ",
    )
