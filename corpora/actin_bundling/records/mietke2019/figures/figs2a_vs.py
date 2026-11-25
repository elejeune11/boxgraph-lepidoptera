import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

def _cal_figs2a(img):
    box = C.find_plot_box(img)
    c0, r0, c1, r1 = box
    x_ticks = C.detect_interior_ticks(img, box, "bottom", band_lo=10, band_hi=20)
    y_ticks = C.detect_interior_ticks(img, box, "left", band_lo=10, band_hi=22)
    xa = C.match_ticks_log10(
        x_ticks,
        [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0],
        10.0, 1000.0, c0, c1, tol_frac=0.01,
    )
    ya = C.match_ticks_log10(y_ticks, [1e-8, 1e-6, 1e-4, 1e-2], 1e-9, 1e-2, r1, r0)
    excl = [(c0, r0 + 1000, c0 + 780, r1)]
    return box, xa, ya, excl

def read(img, at):
    box, xa, ya, excl = _cal_figs2a(img)
    return C.read_markers_xo(img, xa, ya, box, kind="o", exclude_rects=excl)
