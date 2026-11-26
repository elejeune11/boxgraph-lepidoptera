"""Calibrated reader for Fig. 1(d) (p. 2): espin:actin SAXS intensity I(q) at R = 0.05 (blue, upper panel)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import calibrate_fig1d, sample_curve_from_mask


def read(image, at):
    x_axis, y_axis, (cL, rT, cR, rB), lr_diff_px = calibrate_fig1d(image, "upper")
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    # Dark blue/indigo curve (R = 0.05) in upper panel
    blue = (b > 80.0) & (b - r > 30.0) & (b - g > 25.0)
    return sample_curve_from_mask(
        blue, x_axis, y_axis, at, r_lo=rT + 100, r_hi=int(round(y_axis.to_pixel(50.0))),
        extra_err_px=0.5 * lr_diff_px, search_dx=30,
        note_prefix=f"fig1d espin R=0.05 blue (upper, lr_diff={lr_diff_px:.1f}px); ",
    )
