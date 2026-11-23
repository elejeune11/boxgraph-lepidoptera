"""Calibrated reader for Main Fig. 2A bottom-right (p. 3): FRAP experimental variance with error bars."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def calibrate_fig2a_diff(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=140.0)
    h_top, h_bot = hl[0], hl[-1]
    v_box = [v for v in vl if v[0] > 40.0]
    v_left, v_right = v_box[0], v_box[-1]
    _, xt_bot, _, yt_right = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 15), margin=12, max_width=10)
    # 4 interior x-ticks at 20, 40, 60, 80 (ignore pre-bleach red ticks near bottom left)
    xt_cand = [t for t in xt_bot if t > v_left[0] + 0.14 * (v_right[0] - v_left[0])]
    xt_int = match_ticks_regular(xt_cand, 4, tol=3.0)
    x_axis = verify_axis("linear", [(v_left[0], 0.0)] + list(zip(xt_int, [20.0, 40.0, 60.0, 80.0])) + [(v_right[0], 100.0)], max_residual_px=1.5)
    yt_int = match_ticks_regular(yt_right, 8, tol=2.5)
    y_axis = verify_axis("linear", [(h_top[0], 9.0)] + list(zip(yt_int, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])) + [(h_bot[0], 0.0)], max_residual_px=1.5)
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_fig2a_diff(image)
    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180) & (rgb[:, :, 1] < 90) & (rgb[:, :, 2] < 90)
    m = np.zeros_like(red)
    m[int(rt) + 12 : int(rb) - 8, int(cl) + 12 : int(cr) - 6] = red[int(rt) + 12 : int(rb) - 8, int(cl) + 12 : int(cr) - 6]
    # Exclude legend in top right
    m[: int(rt) + 210, int(cl) + int(0.25 * (cr - cl)) :] = False

    # Group red error bars by x-column clusters (merges bars crossed by the green model line)
    x_groups = clusters(np.where(m.any(axis=0))[0], gap=25)
    cap_err_px = 2.0 + y_axis.residual_px
    mid_err_px = 3.0 + y_axis.residual_px

    out = []
    for _, c0, c1 in x_groups:
        rr, cc = np.nonzero(m[:, c0 : c1 + 1])
        if len(rr) < 25 or (rr.max() - rr.min()) < 8:
            continue  # skip tiny pre-bleach ticks on baseline (y=0, no error bar)
        cc = cc + c0
        c_mid = 0.5 * (float(cc.min()) + float(cc.max())) + 0.5
        r_top_cap = float(rr.min()) + 0.5
        r_bot_cap = float(rr.max()) + 0.5
        r_mid = 0.5 * (r_top_cap + r_bot_cap)

        x_val = float(x_axis.to_value(c_mid))
        y_val, lo, hi = value_interval(y_axis, r_mid, mid_err_px)
        b_hi = float(y_axis.to_value(r_top_cap))
        b_lo = float(y_axis.to_value(r_bot_cap))
        b_hi_o = float(y_axis.to_value(r_top_cap - cap_err_px))
        b_lo_o = float(y_axis.to_value(r_bot_cap + cap_err_px))
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "bar_lo": min(b_lo, b_hi),
            "bar_hi": max(b_lo, b_hi),
            "bar_lo_outer": min(b_lo_o, b_hi_o),
            "bar_hi_outer": max(b_lo_o, b_hi_o),
            "px": [c_mid, r_mid],
            "px_lo": r_top_cap,
            "px_hi": r_bot_cap,
            "note": f"error bar caps ({r_top_cap:.1f}, {r_bot_cap:.1f}) px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return sorted(out, key=lambda r: r["x"])
