import importlib.util
from pathlib import Path
import numpy as np
from scipy import ndimage

_spec = importlib.util.spec_from_file_location("_common", Path(__file__).resolve().parent / "_common.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)


def read(img, at):
    box = C.find_plot_box(img)
    c0, r0, c1, r1 = box
    x_ticks = C.detect_interior_ticks(img, box, "bottom")
    y_ticks = C.detect_interior_ticks(img, box, "left")
    xa = C.match_ticks_linear(x_ticks, [0.0, 0.25, 0.5, 0.75, 1.0], 0.0, 1.0, c0, c1)
    ya = C.match_ticks_linear(y_ticks, [0.0, 5.0, 10.0, 15.0, 20.0], 0.0, 20.0, r1, r0)
    mask = C.rgb_mask(img, (0, 114, 189), 65.0)
    H, W = mask.shape
    roi = np.zeros_like(mask)
    roi[int(r0) + 6:min(H, int(r1) + 16), int(c0) + 8:int(c1) - 8] = True
    # Exclude legend box (s/L > 0.56, c in [3, 12])
    roi[int(ya.to_pixel(12.0)):int(ya.to_pixel(3.0)), int(xa.to_pixel(0.56)):] = False
    m = mask & roi
    # Bridge the dark bottom-border stroke (r1-8..r1+8) where the thick blue c->0 curve passes under it
    b_lo = max(0, int(r1) - 24)
    b_hi = min(H, int(r1) + 16)
    m[b_lo:b_hi, :] = ndimage.binary_closing(m[b_lo:b_hi, :], structure=np.ones((22, 1)))
    rows = C.read_vertical_curve(m, xa, ya, at, box)
    # Where the curve touches the top edge of the opaque black bottom border (y ~ 0.1..0.4)
    # without emerging below it, its lower half is occluded behind the border stroke (~16 px);
    # extend the center down by the missing half-thickness and widen the interval to r1.
    border_top = int(r1) - 10
    for r in rows:
        if r.get("y") is None:
            continue
        col = int(round(r["px"][0]))
        idx = np.where(m[:, col])[0]
        if len(idx) and idx[-1] >= border_top and idx[-1] < int(r1) + 2:
            # Curve enters the opaque border stroke from above; true center is between idx[0]+15 and r1
            est_ctr = min(float(r1), float(idx[0]) + 15.0)
            lo_px = float(idx[0])
            hi_px = float(r1) + 2.0
            yv = ya.to_value(est_ctr)
            a, b = ya.to_value(lo_px), ya.to_value(hi_px)
            r["y"] = float(yv)
            r["lo"] = float(min(a, b))
            r["hi"] = float(max(a, b))
            r["px"][1] = est_ctr
            r["px_lo"] = lo_px
            r["px_hi"] = hi_px
            r["note"] = f"partially occluded by bottom border (cal={ya.residual_px:.2f} px)"
    return rows
