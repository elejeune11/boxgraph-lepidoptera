"""Calibrated reader for Fig. 1(d) (p. 2): pure F-actin SAXS intensity I(q) (green, lower panel)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import calibrate_fig1d, sample_curve_from_mask


def read(image, at):
    x_axis, y_axis, (cL, rT, cR, rB), lr_diff_px = calibrate_fig1d(image, "lower")
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    # Green curve (pure F-actin) in lower panel below '//' break
    green = (g > 60.0) & (g - r > 20.0) & (g - b > 15.0)
    return sample_curve_from_mask(
        green, x_axis, y_axis, at, r_lo=int(round(y_axis.to_pixel(46.0))), r_hi=rB - 5,
        extra_err_px=0.5 * lr_diff_px, search_dx=32,
        note_prefix=f"fig1d pure F-actin green (lower, lr_diff={lr_diff_px:.1f}px); ",
    )
