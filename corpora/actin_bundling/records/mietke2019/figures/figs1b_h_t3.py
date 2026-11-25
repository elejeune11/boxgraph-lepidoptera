import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

def read(img, at):
    box = C.find_plot_box(img)
    c0, r0, c1, r1 = box
    xa = C.match_ticks_linear(C.detect_interior_ticks(img, box, "bottom"), [0.0, 0.25, 0.5, 0.75, 1.0], 0.0, 1.0, c0, c1)
    ya = C.match_ticks_linear(C.detect_interior_ticks(img, box, "left"), [0.6, 0.7, 0.8, 0.9, 1.1, 1.2], 0.6, 1.2, r1, r0)
    mask = C.rgb_mask(img, (161, 218, 180), 45.0) & C.box_interior_mask(img.shape[:2], box, inset=8)
    # Exclude text label "t3" below curve near u/L0 in [0.04, 0.15], h < 0.70
    mask[int(ya.to_pixel(0.70)):, int(xa.to_pixel(0.04)):int(xa.to_pixel(0.15))] = False
    return C.read_vertical_curve(mask, xa, ya, at, box)
