"""Calibrated reader for Fig. 6C (p. 9): 4 nM CP, 0 nM VASP fit curve (K_obs = 0.049 s^-1)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_fig6c(image)
    def mask_fn(m, xax, yax):
        m[int(yax.to_pixel(0.51)) :, int(xax.to_pixel(140.0)) :] = False
    return C.trace_monotonic_capping_curve(image, xa, ya, fr, seed_x=25.0, seed_rank=0, at=at, mask_fn=mask_fn, seed_thresh=70.0, march_thresh=80.0)
