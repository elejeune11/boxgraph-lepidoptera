"""Calibrated reader for Fig. 3G (p. 5): 1 - CF curve for p(ADP-Pi) < 0.5 (open circles / upper curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig3g(image)
    mask = C.black_ink_mask(image, gray_max=155.0, sat_max=35.0)
    mask[: int(ya.to_pixel(0.42)), int(xa.to_pixel(1.2)) :] = False
    return C.trace_curve_at(mask, xa, ya, fr["box"], at, half_w=3, inset_px=2, extra_err_px=1.2, select="upper")
