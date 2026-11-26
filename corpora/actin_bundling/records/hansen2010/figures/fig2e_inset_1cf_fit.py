"""Calibrated reader for Fig. 2E inset (p. 4): single-exponential fit of 1 - CF (tau_1 = 0.44 s)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig2e_inset(image)
    mask = C.red_mask(image, r_min=115.0, rg_diff=35.0, rb_diff=35.0)
    return C.trace_curve_at(mask, xa, ya, fr["box"], at, half_w=2, inset_px=2, extra_err_px=1.0)
