"""Calibrated reader for Fig. 2H (p. 4): VASP LIL->AAA (red dash-dot 1 - CF curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig2h(image)
    mask = C.red_mask(image, r_min=110.0, rg_diff=35.0, rb_diff=35.0)
    # Exclude red legend swatch (x > 0.9, y > 0.30)
    mask[: int(ya.to_pixel(0.30)), int(xa.to_pixel(0.90)) :] = False
    return C.trace_curve_at(mask, xa, ya, fr["box"], at, half_w=3, inset_px=2, extra_err_px=1.2, search_w_fallback=12)
