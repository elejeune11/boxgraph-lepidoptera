"""Calibrated reader for Fig. S3E (p. 17): green dashed [ADP-Pi] nucleotide-state probability model curve."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS3e(image)
    mask = C.green_mask(image, g_min=75.0, gr_diff=18.0, gb_diff=15.0)
    mask[: int(ya.to_pixel(0.55)), int(xa.to_pixel(1200.0)) :] = False
    cL, rT, cR, _ = fr["box"]
    return C.trace_curve_at(mask, xa, ya, (cL, rT, cR, fr["x_row"] + 2.0), at, half_w=1, inset_px=0, extra_err_px=1.0, search_w_fallback=14)
