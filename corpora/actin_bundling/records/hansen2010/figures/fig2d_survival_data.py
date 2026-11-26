"""Calibrated reader for Fig. 2D (p. 4): black empirical 1 - CF dwell-time distribution."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig2d(image)
    mask = C.black_ink_mask(image, gray_max=115.0, sat_max=35.0)
    # Exclude annotation text (tau_1, tau_2, R^2, n = 1057) located above curve at x > 1.8 s
    mask[: int(ya.to_pixel(0.16)), int(xa.to_pixel(1.8)) :] = False
    cL, rT, cR, rB = fr["box"]
    rB_eff = min(rB, ya.to_pixel(-0.02))
    return C.trace_curve_at(mask, xa, ya, (cL, rT, cR, rB_eff), at, half_w=3, inset_px=2, extra_err_px=1.2, search_w_fallback=12)
