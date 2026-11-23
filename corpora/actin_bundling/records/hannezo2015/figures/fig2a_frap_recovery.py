"""Calibrated reader for Main Fig. 2A bottom-middle (p. 3): FRAP recovery red '+' data points."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def calibrate_fig2a_tau(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=140.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 15), margin=12, max_width=10)
    xt_int = match_ticks_regular(xt_top, 4, tol=3.0)
    x_axis = verify_axis("linear", [(v_left[0], 0.0)] + list(zip(xt_int, [20.0, 40.0, 60.0, 80.0])) + [(v_right[0], 100.0)], max_residual_px=1.5)
    yt_int = match_ticks_regular(yt_left, 4, tol=3.0)
    y_axis = verify_axis("linear", [(h_top[0], 1.0)] + list(zip(yt_int, [0.8, 0.6, 0.4, 0.2])) + [(h_bot[0], 0.0)], max_residual_px=1.5)
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_fig2a_tau(image)
    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180) & (rgb[:, :, 1] < 90) & (rgb[:, :, 2] < 90)
    m = np.zeros_like(red)
    m[int(rt) + 12 : int(rb) - 4, int(cl) + 8 : int(cr) - 8] = red[int(rt) + 12 : int(rb) - 4, int(cl) + 8 : int(cr) - 8]
    # Exclude legend '+' in upper right
    m[: int(rt) + 190, int(cl) + int(0.35 * (cr - cl)) :] = False

    # Group red '+' pixels by horizontal column clusters (merges '+' markers split by the green curve)
    col_has = np.where(m.any(axis=0))[0]
    x_groups = clusters(col_has, gap=18)

    half_px = 3.0 + y_axis.residual_px
    out = []
    for _, c0, c1 in x_groups:
        rr, cc = np.nonzero(m[:, c0 : c1 + 1])
        if len(rr) < 20:
            continue
        cc = cc + c0
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
