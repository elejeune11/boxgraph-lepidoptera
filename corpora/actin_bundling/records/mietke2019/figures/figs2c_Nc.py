import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

def _cal_figs2c(img):
    box = C.find_plot_box(img)
    c0, r0, c1, r1 = box
    x_ticks = C.detect_interior_ticks(img, box, "bottom", band_lo=10, band_hi=22)
    y_ticks = C.detect_interior_ticks(img, box, "right", band_lo=10, band_hi=22)
    xa = C.match_ticks_linear(x_ticks, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], 0.0, 1.0, c0, c1)
    ya = C.match_ticks_log10(y_ticks, [1e-8, 1e-7, 1e-6, 1e-5, 1e-4], 1e-8, 1e-4, r1, r0)
    return box, xa, ya

def read(img, at):
    box, xa, ya = _cal_figs2c(img)
    # Dashed black molecule-number-error curve lies below y = 1.5e-6
    mask = C.neutral_dark_mask(img, gray_max=100.0, sat_max=30.0) & C.box_interior_mask(img.shape[:2], box, inset=32)
    mask[:int(ya.to_pixel(1.5e-6)), :] = False
    # Exclude legend in bottom-right (t/tau > 0.66, y < 1.8e-7)
    mask[int(ya.to_pixel(1.8e-7)):, int(xa.to_pixel(0.66)):] = False
    return C.read_vertical_curve(mask, xa, ya, at, box, dashed=True, dash_window=55, min_run_px=2)
