"""Calibrated reader for Fig. S3E (p. 17): blue dotted [ADP] nucleotide-state probability model curve."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS3e(image)
    mask = C.blue_mask(image, b_min=95.0, br_diff=25.0, bg_diff=10.0)
    mask[int(ya.to_pixel(0.85)) : int(ya.to_pixel(0.55)), int(xa.to_pixel(1200.0)) :] = False
    cL, rT, cR, _ = fr["box"]
    return C.trace_curve_at(mask, xa, ya, (cL, rT, cR, fr["x_row"] + 2.0), at, half_w=4, inset_px=0, extra_err_px=1.2, search_w_fallback=16)
