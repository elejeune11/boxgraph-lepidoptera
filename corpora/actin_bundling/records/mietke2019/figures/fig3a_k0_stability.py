import importlib.util
from pathlib import Path
import numpy as np

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)
from scipy import ndimage

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
    blue = C.rgb_mask(img, (0, 99, 163), 60.0) & C.box_interior_mask(img.shape[:2], box, inset=6)
    # Exclude legend text "kappa_tilde = 0" (L in [0.52, 0.94], alpha < 1.8)
    blue[int(a_ax.to_pixel(1.8)):, int(L_ax.to_pixel(0.52)):int(L_ax.to_pixel(0.94))] = False
    # Exclude parameter markers (blue cross at 0.82,3.0 and blue circle at 0.95,3.0)
    blue[int(a_ax.to_pixel(3.18)):int(a_ax.to_pixel(2.82)), int(L_ax.to_pixel(0.88)):int(L_ax.to_pixel(1.02))] = False
    # Close thin horizontal gaps where the red dashed curve overlays the blue boundary at alpha > 4.05
    blue = ndimage.binary_closing(blue, structure=np.ones((1, 25)))
    # Select the leftmost boundary run at each alpha
    return C.read_horizontal_curve(blue, a_ax, L_ax, at, box, branch_selector=lambda av, runs, ax: min(runs, key=lambda r: r[0]))
