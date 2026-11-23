"""Calibrated reader for Main Fig. 4F (p. 5): normalized actin intensity vs tube diameter for WT, Cop-/+, and Cop-/- with error bars."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, verify_axis, value_interval


def read(image, at):
    rgb = image.astype(float)
    gray = rgb.mean(axis=2)
    dark = gray < 120.0
    H, W = dark.shape

    # Locate L-shaped axes (ignore top micrograph border rows < 100 and bottom border rows > H-80)
    h_lines = [h for h in clusters(np.where(dark[100 : H - 80].mean(axis=1) > 0.40)[0] + 100)]
    v_lines = [v for v in clusters(np.where(dark[200 : H - 150, : int(0.4 * W)].mean(axis=0) > 0.50)[0])]
    if not h_lines or not v_lines:
        raise RuntimeError("Fig. 4F: could not locate L-shaped axes")
    r_bot, rb0, rb1 = h_lines[-1]
    c_left, cl0, cl1 = v_lines[-1]

    # Leftward y-ticks at 1.5, 1.0, 0.5, 0.0
    yt = [c[0] for c in clusters(np.where(dark[200 : rb1 + 5, cl0 - 28 : cl0 - 6].mean(axis=1) > 0.5)[0] + 200) if (c[2] - c[1] + 1) <= 15]
    if len(yt) != 4:
        raise RuntimeError(f"Fig. 4F: expected 4 y-ticks (1.5, 1.0, 0.5, 0.0), found {len(yt)}: {yt}")
    y_axis = verify_axis("linear", list(zip(yt, [1.5, 1.0, 0.5, 0.0])), max_residual_px=1.5)

    # Downward x-ticks at 1, 2, 3, 4, 5, 6
    xt = [c[0] for c in clusters(np.where(dark[rb1 + 6 : rb1 + 26, cl0 - 5 : int(0.82 * W)].mean(axis=0) > 0.5)[0] + cl0 - 5) if (c[2] - c[1] + 1) <= 15]
    if len(xt) != 6:
        raise RuntimeError(f"Fig. 4F: expected 6 x-ticks (1..6), found {len(xt)}: {xt}")
    x_axis = verify_axis("linear", list(zip(xt, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])), max_residual_px=1.5)

    # Isolate plot interior (strictly left of legend at x ~ 5.8)
    c_max_data = int(x_axis.to_pixel(5.5))
    nw = np.zeros((H, W), bool)
    r0, r1 = int(yt[0]) - 20, int(rb0) - 10
    c0, c1 = int(cl1) + 12, c_max_data
    nw[r0:r1, c0:c1] = (rgb[r0:r1, c0:c1, 0] < 230) | (rgb[r0:r1, c0:c1, 1] < 230) | (rgb[r0:r1, c0:c1, 2] < 230)

    lab, n = ndimage.label(ndimage.binary_dilation(nw, iterations=3))
    cap_err_px = 2.0 + y_axis.residual_px
    mid_err_px = 2.5 + y_axis.residual_px

    out = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero((lab == i) & nw)
        if len(rr) < 500:
            continue
        # Crop tightly around component to run fast 20x20 opening for the filled square marker
        ra, rb_c = int(rr.min()), int(rr.max()) + 1
        ca, cb_c = int(cc.min()), int(cc.max()) + 1
        sub = ((lab[ra:rb_c, ca:cb_c] == i) & nw[ra:rb_c, ca:cb_c])
        sq = ndimage.binary_opening(sub, structure=np.ones((20, 20), bool))
        srr, scc = np.nonzero(sq)
        if len(srr) < 200:
            continue
        r_c = float(srr.mean()) + ra + 0.5
        c_c = float(scc.mean()) + ca + 0.5

        r_top_cap = float(rr.min()) + 0.5
        r_bot_cap = float(rr.max()) + 0.5
        c_left_cap = float(cc.min()) + 0.5
        c_right_cap = float(cc.max()) + 0.5

        x_val = float(x_axis.to_value(c_c))
        x_lo = float(x_axis.to_value(c_left_cap))
        x_hi = float(x_axis.to_value(c_right_cap))
        y_val, lo, hi = value_interval(y_axis, r_c, mid_err_px)
        b_hi = float(y_axis.to_value(r_top_cap))
        b_lo = float(y_axis.to_value(r_bot_cap))
        b_hi_o = float(y_axis.to_value(r_top_cap - cap_err_px))
        b_lo_o = float(y_axis.to_value(r_bot_cap + cap_err_px))

        rgb_med = [int(np.median(rgb[srr + ra, scc + ca, k])) for k in range(3)]
        if max(rgb_med) < 50:
            lbl = "Wild-type"
        elif rgb_med[0] < 80:
            lbl = "Cop-/+"
        else:
            lbl = "Cop-/-"

        out.append({
            "x": x_val,
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "bar_lo": min(b_lo, b_hi),
            "bar_hi": max(b_lo, b_hi),
            "bar_lo_outer": min(b_lo_o, b_hi_o),
            "bar_hi_outer": max(b_lo_o, b_hi_o),
            "px": [c_c, r_c],
            "px_lo": r_top_cap,
            "px_hi": r_bot_cap,
            "note": f"{lbl}: x_bar=[{x_lo:.3f}, {x_hi:.3f}] um (+/-{0.5*(x_hi-x_lo):.3f} um); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return sorted(out, key=lambda r: r["x"])
