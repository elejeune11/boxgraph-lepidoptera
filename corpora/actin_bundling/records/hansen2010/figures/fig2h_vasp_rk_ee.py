"""Calibrated reader for Fig. 2H (p. 4): VASP RK->EE (blue dotted 1 - CF curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig2h(image)
    mask = C.blue_mask(image, b_min=90.0, br_diff=20.0, bg_diff=10.0)
    # Exclude blue legend swatch (x > 0.9, y > 0.25)
    mask[: int(ya.to_pixel(0.25)), int(xa.to_pixel(0.90)) :] = False
    return C.trace_curve_at(mask, xa, ya, fr["box"], at, half_w=4, inset_px=2, extra_err_px=1.5, search_w_fallback=14)
