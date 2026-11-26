"""Calibrated reader for Fig. 3I (p. 5): red dashed p(ADP-Pi) nucleotide-state model curve."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya_r, fr = C.calibrate_fig3i(image)
    mask = C.red_mask(image, r_min=110.0, rg_diff=35.0, rb_diff=35.0)
    # Exclude red legend entry (y > 0.50, x > 8)
    mask[: int(ya_r.to_pixel(0.50)), int(xa.to_pixel(8.0)) :] = False
    return C.trace_curve_at(mask, xa, ya_r, fr["box"], at, half_w=5, inset_px=3, extra_err_px=1.2, search_w_fallback=14)
