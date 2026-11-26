"""Calibrated reader for Fig. 4E inset (p. 7): single-exponential fit of 1 - CF (tau_1 = 1.45 s)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig4e_inset(image)
    mask = C.red_mask(image, r_min=95.0, rg_diff=25.0, rb_diff=25.0)
    cL, rT, cR, _ = fr["box"]
    return C.trace_curve_at(mask, xa, ya, (cL, rT, cR, fr["x_row"] + 2.0), at, half_w=4, inset_px=0, extra_err_px=1.0, search_w_fallback=22)
