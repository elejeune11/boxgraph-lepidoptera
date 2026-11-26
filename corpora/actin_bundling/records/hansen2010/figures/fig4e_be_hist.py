"""Calibrated reader for Fig. 4E (p. 7): barbed-end dwell-time histogram (0.2 s red bars)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from boxgraph import digitize as D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    fr = C.find_plot_frame(image)
    yt = C.detect_spine_ticks(image, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = C.match_linear_axis(yt, [0, 50, 100, 150, 200], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 14.0)])
    rm = C.red_mask(image, r_min=135.0, rg_diff=45.0, rb_diff=45.0)
    # Mask out the inset plot box in upper right (x > 2.4 s, y > 40)
    rm[: int(ya.to_pixel(40.0)), int(xa.to_pixel(2.4)) :] = False
    rT, rB = int(round(fr["rT"])), int(round(fr["rB"]))
    cal_err = max(ya.residual_px, 0.5)
    out = []
    # Read 0.2 s bins from 0.1 to 7.9 s (40 bins covering almost all n = 1166 events)
    for k in range(40):
        x_lo, x_hi = k * 0.2, (k + 1) * 0.2
        x_mid = 0.5 * (x_lo + x_hi)
        c0 = int(round(xa.to_pixel(x_lo))) + 3
        c1 = int(round(xa.to_pixel(x_hi))) - 3
        c_mid = xa.to_pixel(x_mid)
        sub = rm[rT + 2 : rB - 1, c0 : max(c0 + 1, c1)]
        ys = np.nonzero(sub.any(axis=1))[0]
        if len(ys) == 0:
            r_top = fr["x_row"]
            _, _, hi = C.val_interval(ya, r_top, 1.5 + cal_err)
            out.append({
                "x": round(x_mid, 3), "y": 0.0, "lo": 0.0, "hi": hi,
                "px": [float(c_mid), float(r_top)],
                "px_lo": float(r_top - 2.0), "px_hi": float(r_top),
                "note": f"bin [{x_lo:.1f},{x_hi:.1f}]s empty",
            })
        else:
            col_tops = [rT + 2 + int(np.nonzero(sub[:, j])[0][0]) for j in range(sub.shape[1]) if sub[:, j].any()]
            r_top = float(np.median(col_tops))
            y, lo, hi = C.val_interval(ya, r_top, 1.5 + cal_err)
            lo = max(0.0, lo)
            out.append({
                "x": round(x_mid, 3), "y": y, "lo": lo, "hi": hi,
                "px": [float(c_mid), float(r_top)],
                "px_lo": float(r_top - 2.0), "px_hi": float(r_top + 2.0),
                "note": f"bin [{x_lo:.1f},{x_hi:.1f}]s red bar",
            })
    return out
