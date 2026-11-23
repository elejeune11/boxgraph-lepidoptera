"""Calibrated reader for SI Fig. S9 (p. 16): zero-friction bifurcation onset abscissae chi_i^(0)/phi for modes i = 1..4."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def calibrate_figs9(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=130.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    _, xt_bot, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(6, 22), margin=15, max_width=14)

    # X ticks: left = 0.5, 6 interior ticks = 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, right = 4.0
    xt_int = match_ticks_regular(xt_bot, 6, tol=3.5)
    x_marks = [(v_left[0], 0.5)] + list(zip(xt_int, [1.0, 1.5, 2.0, 2.5, 3.0, 3.5])) + [(v_right[0], 4.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=3.0)

    # Y ticks: top = 7, 6 interior ticks = 6, 5, 4, 3, 2, 1, bottom = 0
    yt_int = match_ticks_regular(yt_left, 6, tol=3.5)
    y_marks = [(h_top[0], 7.0)] + list(zip(yt_int, [6.0, 5.0, 4.0, 3.0, 2.0, 1.0])) + [(h_bot[0], 0.0)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=3.0)
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs9(image)
    rgb = image.astype(float)
    blue = (rgb[:, :, 2] > 180.0) & (rgb[:, :, 0] < 65.0) & (rgb[:, :, 1] < 65.0)

    # Measure each branch just above the homogeneous line y = 1.0 (in band y in [1.05, 1.25])
    # and extrapolate linearly in (y - 1)^2 -> 0 (pitchfork: chi/phi = chi0/phi + c*(y - 1) near y = 1)
    r_one = float(y_axis.to_pixel(1.0))
    r_lo = int(round(y_axis.to_pixel(1.30)))
    r_hi = int(round(y_axis.to_pixel(1.05)))

    # Identify 4 distinct branch x-bands at y = 1.15
    r_ref = int(round(y_axis.to_pixel(1.15)))
    ref_cols = clusters(np.where(blue[r_ref - 15 : r_ref + 16, int(cl) + 20 : int(cr) - 20].any(axis=0))[0] + int(cl) + 20, gap=30)
    if len(ref_cols) != 4:
        raise RuntimeError(f"Fig. S9: expected 4 bifurcating branches near y=1.15, found {len(ref_cols)}: {ref_cols}")

    out = []
    for mode_i, (c_ref, _, _) in enumerate(ref_cols, start=1):
        # Collect blue pixels of this branch for y in [1.04, 1.35] (excluding horizontal y=1 dashed line!)
        c_win0 = max(int(cl) + 10, int(c_ref) - 120)
        c_win1 = min(int(cr) - 10, int(c_ref) + 120)
        ys_v = []
        xs_c = []
        for r in range(r_lo, r_hi + 1):
            cs = np.where(blue[r, c_win0:c_win1])[0] + c_win0
            if len(cs) > 0:
                ys_v.append(float(y_axis.to_value(r + 0.5)))
                xs_c.append(float(cs.mean()) + 0.5)
        ys_v = np.array(ys_v)
        xs_c = np.array(xs_c)
        # Fit x_px(y) linearly over y in [1.05, 1.30] and extrapolate to y = 1.0
        slope, intercept = np.polyfit(ys_v - 1.0, xs_c, 1)
        c_onset = float(intercept)
        resid_fit = float(np.std(xs_c - (slope * (ys_v - 1.0) + intercept)))
        half_px = max(4.0, 2.0 * resid_fit) + x_axis.residual_px + 2.0
        val, lo, hi = value_interval(x_axis, c_onset, half_px)
        out.append({
            "x": float(mode_i),
            "y": val,
            "lo": lo,
            "hi": hi,
            "px": [c_onset, r_one],
            "note": f"mode i={mode_i} onset extrapolated from y in [1.05,1.30] ({len(xs_c)} rows); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return out
