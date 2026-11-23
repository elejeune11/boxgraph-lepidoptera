"""Calibrated reader for SI Fig. S15 (p. 22): steady-state number of rings vs chi (blue markers with error bars)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def calibrate_figs15(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=130.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(6, 22), margin=15, max_width=14)

    # X ticks: left border = 0, 5 interior ticks on top = 1, 2, 3, 4, 5
    xt_int = match_ticks_regular(xt_top, 5, tol=3.5)
    x_axis = verify_axis("linear", [(v_left[0], 0.0)] + list(zip(xt_int, [1.0, 2.0, 3.0, 4.0, 5.0])), max_residual_px=2.5)

    # Y ticks: 7 interior ticks = 35, 30, 25, 20, 15, 10, 5 (exclude marker-thickened bottom border)
    yt_int = match_ticks_regular(yt_left, 7, tol=3.5)
    y_axis = verify_axis(
        "linear",
        list(zip(yt_int, [35.0, 30.0, 25.0, 20.0, 15.0, 10.0, 5.0])),
        max_residual_px=1.5,
    )
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs15(image)
    rgb = image.astype(float)
    # Blue markers & error bars (note red curve crosses near the centers of the blue error bars)
    blue = (rgb[:, :, 2] > 160.0) & (rgb[:, :, 0] < 95.0) & (rgb[:, :, 1] < 95.0)
    r_lo, r_hi = int(rt) + 6, int(rb) - 25
    c_lo, c_hi = int(cl) + 20, int(cr) - 10

    # Locate the 40 vertical error-bar stems by counting blue pixels per column (stems have >= 35 blue px vertically,
    # whereas horizontal caps only contribute ~4-6 px in off-stem columns)
    col_counts = blue[r_lo:r_hi, c_lo:c_hi].sum(axis=0)
    stem_cols = clusters(np.where(col_counts >= 35)[0] + c_lo, gap=12)

    cap_err_px = 2.5 + y_axis.residual_px
    mid_err_px = 3.0 + y_axis.residual_px
    out = []
    for c_mid, c0, c1 in stem_cols:
        # Look within +/- 8 px of stem center for the top and bottom cap rows
        cc0 = max(c_lo, int(round(c_mid)) - 8)
        cc1 = min(c_hi, int(round(c_mid)) + 9)
        rr = np.where(blue[r_lo:r_hi, cc0:cc1].any(axis=1))[0] + r_lo
        if len(rr) < 20:
            continue
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
            "note": f"blue stem at col {c_mid:.1f}, caps ({r_top_cap:.1f},{r_bot_cap:.1f})px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return sorted(out, key=lambda r: r["x"])
