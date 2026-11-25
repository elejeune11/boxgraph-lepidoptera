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
    y_ticks = C.detect_interior_ticks(img, box, "left")
    xa = C.match_ticks_linear(x_ticks, [0.0, 0.25, 0.5, 0.75, 1.0], 0.0, 1.0, c0, c1)
    ya = C.match_ticks_linear(y_ticks, [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0], -1.0, 5.0, r1, r0)
    mask = C.rgb_mask(img, (0, 114, 189), 65.0) & C.box_interior_mask(img.shape[:2], box, inset=8)
    # Exclude legend in top-left (s/L < 0.45, y > 1.8)
    mask[:int(ya.to_pixel(1.8)), :int(xa.to_pixel(0.45))] = False
    # Pick the thick solid curve rather than the thin horizontal dashed line at y = 0
    def pick(xv, runs, ax):
        thick = [r for r in runs if 2.0 * r[1] >= 18.0]
        if thick:
            return max(thick, key=lambda r: r[1])
        return max(runs, key=lambda r: abs(ax.to_value(r[0])))
    return C.read_vertical_curve(mask, xa, ya, at, box, branch_selector=pick)
