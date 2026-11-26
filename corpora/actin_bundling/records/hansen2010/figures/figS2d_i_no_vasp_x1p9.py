"""Calibrated reader for Fig. S2D(i) (p. 16): -VASP (x 1.9) dashed curve at 100 mM KCl."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    xa, ya, fr = C.calibrate_figS2d_i(image)
    m = C.gray(image) < 160.0
    # Exclude upper-right legend and '[KCl]=100mM' label
    m[: int(ya.to_pixel(4000.0)), int(xa.to_pixel(5.0)) :] = False
    m[: int(ya.to_pixel(4900.0)), :] = False
    def pick(x, cl):
        if len(cl) < 3:
            return None
        # At x=1, 2, 4, the dashed line sits strictly between +50 nM VASP (top) and -VASP (bottom)
        # If multiple dashes appear inside half_w=10 at steep x=1, pick the middle cluster of the interior
        interior = list(range(1, len(cl) - 1))
        return interior[len(interior) // 2]
    return C.read_ranked_dark_curves(
        image, xa, ya, fr, at, pick,
        mask=m, half_w=10, merge_gap=6, y_min_val=50.0, y_max_val=4800.0,
    )
