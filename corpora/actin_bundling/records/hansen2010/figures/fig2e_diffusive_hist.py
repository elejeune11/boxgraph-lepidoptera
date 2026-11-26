"""Calibrated reader for Fig. 2E (p. 4): diffusive dwell-time histogram (0.1 s red bars)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from boxgraph import digitize as D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    fr = C.find_plot_frame(image)
    # Fig. 2E histogram frame spans x = [0, 5] s and y = [0, 250] directly (flush frame, 50 bins of 0.1 s)
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 5.0)])
    ya = D.Axis("linear", [(fr["x_row"], 0.0), (fr["top_row"], 250.0)])
    rm = C.red_mask(image, r_min=135.0, rg_diff=45.0, rb_diff=45.0)
    rT, rB = int(round(fr["rT"])), int(round(fr["rB"]))
    # Mask out the inset plot box in upper right (x > 1.3 s, y > 45)
    rm[: int(ya.to_pixel(45.0)), int(xa.to_pixel(1.3)) :] = False
    out = []
    # Read 0.1 s bins from [0.0, 0.1] up to [2.9, 3.0] (center 0.05..2.95)
    for k in range(30):
        x_lo, x_hi = k * 0.1, (k + 1) * 0.1
        x_mid = 0.5 * (x_lo + x_hi)
        c0 = int(round(xa.to_pixel(x_lo))) + 3
        c1 = int(round(xa.to_pixel(x_hi))) - 3
        c_mid = xa.to_pixel(x_mid)
        sub = rm[rT + 2 : rB - 1, c0 : max(c0 + 1, c1)]
        ys = np.nonzero(sub.any(axis=1))[0]
        if len(ys) == 0:
            r_top = fr["x_row"]
            y, lo, hi = C.val_interval(ya, r_top, 2.0)
            lo = max(0.0, lo)
            out.append({
                "x": round(x_mid, 3), "y": 0.0, "lo": 0.0, "hi": hi,
                "px": [float(c_mid), float(r_top)],
                "px_lo": float(r_top - 2.0), "px_hi": float(r_top),
                "note": f"bin [{x_lo:.1f},{x_hi:.1f}]s empty",
            })
        else:
            # Top edge of the red bar (use median across columns to avoid 1-px JPEG fringe)
            col_tops = [rT + 2 + int(np.nonzero(sub[:, j])[0][0]) for j in range(sub.shape[1]) if sub[:, j].any()]
            r_top = float(np.median(col_tops))
            y, lo, hi = C.val_interval(ya, r_top, 2.0)
            lo = max(0.0, lo)
            out.append({
                "x": round(x_mid, 3), "y": y, "lo": lo, "hi": hi,
                "px": [float(c_mid), float(r_top)],
                "px_lo": float(r_top - 2.0), "px_hi": float(r_top + 2.0),
                "note": f"bin [{x_lo:.1f},{x_hi:.1f}]s red bar",
            })
    return out
