import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

def read(img, at):
    box = C.find_plot_box(img)
    c0, r0, c1, r1 = box
    x_ticks = C.detect_interior_ticks(img, box, "bottom")
    r_ticks = C.detect_interior_ticks(img, box, "right")
    xa = C.match_ticks_linear(x_ticks, [0.0, 0.25, 0.5, 0.75, 1.0], 0.0, 1.0, c0, c1)
    ya = C.calibrate_shifted_linear_ticks(r_ticks, [0.5, 1.0, 1.5], ascending_pixel_means_ascending_value=False)
    mask = C.rgb_mask(img, (217, 83, 25), 65.0) & C.box_interior_mask(img.shape[:2], box, inset=8)
    mask[:int(ya.to_pixel(1.3)), int(xa.to_pixel(0.65)):] = False
    return C.read_vertical_curve(mask, xa, ya, at, box)
