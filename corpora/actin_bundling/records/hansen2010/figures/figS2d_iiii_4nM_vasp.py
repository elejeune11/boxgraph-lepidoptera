"""Calibrated reader for Fig. S2D(iiii) (p. 16): 4 nM CP, +50 nM VASP, 100 mM KCl."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS2d_iii(image)
    def mask_fn(m, xax, yax):
        m[int(yax.to_pixel(0.20)) :, int(xax.to_pixel(270.0)) :] = False
    return C.trace_monotonic_capping_curve(image, xa, ya, fr, seed_x=260.0, seed_rank=1, at=at, mask_fn=mask_fn, seed_thresh=70.0, march_thresh=95.0)
