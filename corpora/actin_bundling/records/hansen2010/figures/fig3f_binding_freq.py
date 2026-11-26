"""Calibrated reader for Fig. 3F (p. 5): Cy3-VASP binding frequency (# molecules) vs filament age."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya_l, _, fr = C.calibrate_fig3f(image)
    mask = C.black_ink_mask(image, gray_max=120.0, sat_max=35.0)
    # Exclude legend in upper-left/center (y > 19.5, x < 14)
    mask[: int(ya_l.to_pixel(19.5)), : int(xa.to_pixel(14.5))] = False
    return C.trace_curve_at(mask, xa, ya_l, fr["box"], at, half_w=4, inset_px=3, extra_err_px=1.5)
