import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

def _cal_fig3a(img):
    box = C.find_plot_box(img)
    c0, r0, c1, r1 = box
    x_ticks = C.detect_interior_ticks(img, box, "bottom")
    y_ticks = C.detect_interior_ticks(img, box, "left")
    L_ax = C.match_ticks_linear(x_ticks, [0.0, 0.25, 0.5, 0.75, 1.0, 1.25], 0.0, 1.25, c0, c1)
    a_ax = C.match_ticks_linear(y_ticks, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], 0.0, 5.0, r1, r0)
    return box, L_ax, a_ax

def read(img, at):
    box, L_ax, a_ax = _cal_fig3a(img)
    # Dark blue oscillatory region (alone or overlapped with pink oscillatory region)
    reg = (C.rgb_mask(img, (91, 171, 224), 40.0) | C.rgb_mask(img, (123, 138, 181), 40.0)) & C.box_interior_mask(img.shape[:2], box, inset=6)
    return C.read_horizontal_region_right_edge(reg, a_ax, L_ax, at, box)
