"""Calibrated reader for Main Fig. 3C top-right (p. 4): kkv mutant tube diameter vs actin concentration (scatter '+' points)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def calibrate_fig3c_diam(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=140.0)
    # Exclude title text lines above row 140
    h_box = [h for h in hl if h[0] > 140.0]
    h_top, h_bot = h_box[0], h_box[-1]
    v_left, v_right = vl[0], vl[-1]

    _, xt_bot, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 15), margin=12, max_width=10)
    # X ticks: left = 0, interior = 50, 100, 150, 200, right = 250
    xt_int = match_ticks_regular(xt_bot, 4, tol=3.0)
    x_axis = verify_axis("linear", [(v_left[0], 0.0)] + list(zip(xt_int, [50.0, 100.0, 150.0, 200.0])) + [(v_right[0], 250.0)], max_residual_px=1.5)

    # Y ticks: top = 8, 7 interior ticks = 7, 6, 5, 4, 3, 2, 1, bottom = 0
    yt_int = match_ticks_regular([t for t in yt_left if t > h_top[0] + 20], 7, tol=2.5)
    y_axis = verify_axis("linear", [(h_top[0], 8.0)] + list(zip(yt_int, [7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])) + [(h_bot[0], 0.0)], max_residual_px=1.5)
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_fig3c_diam(image)
    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 170) & (rgb[:, :, 1] < 85) & (rgb[:, :, 2] < 85)
    m = np.zeros_like(red)
    m[int(rt) + 10 : int(rb) - 8, int(cl) + 10 : int(cr) - 10] = red[int(rt) + 10 : int(rb) - 8, int(cl) + 10 : int(cr) - 10]

    # Dilate lightly to rejoin '+' symbols split by the blue trend curve
    lab, n = ndimage.label(ndimage.binary_dilation(m, iterations=2))
    half_px = 3.0 + y_axis.residual_px
    out = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero((lab == i) & m)
        if len(rr) < 12:
            continue
        c_mid = 0.5 * (float(cc.min()) + float(cc.max())) + 0.5
        r_mid = 0.5 * (float(rr.min()) + float(rr.max())) + 0.5
        x_val = float(x_axis.to_value(c_mid))
        y_val, lo, hi = value_interval(y_axis, r_mid, half_px)
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [c_mid, r_mid],
            "px_lo": float(r_mid - half_px),
            "px_hi": float(r_mid + half_px),
            "note": f"red '+' marker ({len(rr)} px); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return sorted(out, key=lambda r: r["x"])
