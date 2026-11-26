"""Calibrated reader for Fig. 2H (p. 4): VASP + 2 uM Actin (green dashed 1 - CF curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig2h(image)
    mask = C.green_mask(image, g_min=65.0, gr_diff=15.0, gb_diff=15.0)
    # Exclude green legend swatch (x > 0.9, y > 0.20)
    mask[: int(ya.to_pixel(0.20)), int(xa.to_pixel(0.90)) :] = False
    return C.trace_curve_at(mask, xa, ya, fr["box"], at, half_w=4, inset_px=2, extra_err_px=1.2, search_w_fallback=14)
