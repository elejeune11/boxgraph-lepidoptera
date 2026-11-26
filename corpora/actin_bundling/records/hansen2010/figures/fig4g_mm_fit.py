"""Calibrated reader for Fig. 4G (p. 7): red offset Michaelis-Menten fit curve (Kd = 9.2 nM)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig4g(image)
    mask = C.red_mask(image, r_min=115.0, rg_diff=35.0, rb_diff=35.0)
    return C.trace_curve_at(mask, xa, ya, fr["box"], at, half_w=2, inset_px=2, extra_err_px=1.0, search_w_fallback=10)
