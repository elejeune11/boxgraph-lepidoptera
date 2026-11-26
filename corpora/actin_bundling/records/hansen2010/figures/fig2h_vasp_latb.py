"""Calibrated reader for Fig. 2H (p. 4): VASP +LatB (black solid 1 - CF curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig2h(image)
    mask = C.black_ink_mask(image, gray_max=115.0, sat_max=35.0)
    # Exclude legend and R^2 label in upper-right quadrant (x > 0.85, y > 0.22)
    mask[: int(ya.to_pixel(0.22)), int(xa.to_pixel(0.85)) :] = False
    cL, rT, cR, rB = fr["box"]
    # Stop above inward x-ticks below y = 0
    rB_eff = min(rB, ya.to_pixel(-0.02))
    return C.trace_curve_at(mask, xa, ya, (cL, rT, cR, rB_eff), at, half_w=2, inset_px=2, extra_err_px=1.0)
