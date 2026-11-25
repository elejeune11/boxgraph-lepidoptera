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
    red = C.rgb_mask(img, (242, 45, 0), 75.0) & C.box_interior_mask(img.shape[:2], box, inset=6)
    # Exclude legend text "kappa_tilde = 0.25" (L in [0.52, 0.98], alpha < 1.6)
    red[int(a_ax.to_pixel(1.6)):, int(L_ax.to_pixel(0.52)):int(L_ax.to_pixel(0.98))] = False
    # Exclude red circle marker at (1.10, 3.0)
    red[int(a_ax.to_pixel(3.22)):int(a_ax.to_pixel(2.78)), int(L_ax.to_pixel(1.03)):int(L_ax.to_pixel(1.17))] = False
    return C.read_horizontal_curve(red, a_ax, L_ax, at, box, dashed=True, dash_window=50)
