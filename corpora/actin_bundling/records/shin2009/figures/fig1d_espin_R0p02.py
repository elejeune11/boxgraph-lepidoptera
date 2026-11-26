"""Calibrated reader for Fig. 1(d) (p. 2): espin:actin SAXS intensity I(q) at R = 0.02 (brown, lower panel)."""
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
    # Brown curve (R = 0.02) in lower panel below '//' break
    brown = (r > 60.0) & (r - b > 30.0) & (r > g + 10.0) & (g > b + 5.0) & (r < 170.0)
    return sample_curve_from_mask(
        brown, x_axis, y_axis, at, r_lo=int(round(y_axis.to_pixel(49.0))), r_hi=rB - 10,
        extra_err_px=0.5 * lr_diff_px, search_dx=32,
        note_prefix=f"fig1d espin R=0.02 brown (lower, lr_diff={lr_diff_px:.1f}px); ",
    )
