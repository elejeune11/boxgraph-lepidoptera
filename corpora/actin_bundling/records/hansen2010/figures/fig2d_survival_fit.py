"""Calibrated reader for Fig. 2D (p. 4): red double-exponential fit curve of 1 - CF."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig2d(image)
    mask = C.red_mask(image, r_min=110.0, rg_diff=35.0, rb_diff=35.0)
    cL, rT, cR, rB = fr["box"]
    rB_eff = min(rB, ya.to_pixel(-0.02))
    return C.trace_curve_at(mask, xa, ya, (cL, rT, cR, rB_eff), at, half_w=2, inset_px=2, extra_err_px=1.0, search_w_fallback=14)
