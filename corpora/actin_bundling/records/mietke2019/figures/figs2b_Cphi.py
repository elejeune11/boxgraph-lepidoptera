import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

import math
def _cal_figs2b(img):
    box = C.find_plot_box(img)
    c0, r0, c1, r1 = box
    # Long major ticks on bottom wall isolate 10^4 and 10^5
    maj_x = C.detect_interior_ticks(img, box, "bottom", band_lo=16, band_hi=28)
    if len(maj_x) < 2:
        raise RuntimeError(f"Expected 2 major x ticks on Fig. S2B, got {maj_x}")
    p1e4, p1e5 = sorted(maj_x)[0], sorted(maj_x)[-1]
    step_dec = p1e5 - p1e4
    # Refine with all top-wall major + minor ticks (8e3..1e5)
    top_x = C.detect_interior_ticks(img, box, "top", band_lo=8, band_hi=18)
    x_marks = [(p1e4, 1e4), (p1e5, 1e5)]
    for val in [8e3, 9e3, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4]:
        pred = p1e4 + step_dec * math.log10(val / 1e4)
        best = min(top_x, key=lambda p: abs(p - pred))
        if abs(best - pred) < 8.0:
            x_marks.append((best, val))
    xa = C.Axis("log10", x_marks)
    if xa.residual_px > 1.5:
        raise RuntimeError(f"Fig. S2B x calibration residual {xa.residual_px:.2f} px")
    # Y axis: interior major ticks at 10^-4, 10^-5, 10^-6 plus bottom border 10^-7
    y_ticks = C.detect_interior_ticks(img, box, "left", band_lo=10, band_hi=22)
    y_sorted = sorted(y_ticks)
    if len(y_sorted) < 3:
        raise RuntimeError(f"Expected 3 interior y ticks on Fig. S2B, got {y_sorted}")
    ya = C.Axis("log10", [(y_sorted[0], 1e-4), (y_sorted[1], 1e-5), (y_sorted[2], 1e-6), (r1, 1e-7)])
    if ya.residual_px > 1.5:
        raise RuntimeError(f"Fig. S2B y calibration residual {ya.residual_px:.2f} px")
    excl = [(c0, r0 + 900, c0 + 700, r1), (c0 + 600, r0, c1, r0 + 280)]
    return box, xa, ya, excl

def read(img, at):
    box, xa, ya, excl = _cal_figs2b(img)
    return C.read_markers_xo(img, xa, ya, box, kind="o", exclude_rects=excl)
