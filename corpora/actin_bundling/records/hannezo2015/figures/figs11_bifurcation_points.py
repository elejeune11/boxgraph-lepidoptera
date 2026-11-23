"""Calibrated reader for SI Fig. S11 (p. 19): primary linear bifurcation thresholds chi_i^(0) along max rho = 1 (alpha = 35, phi = 0.1)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, find_frame, find_inward_ticks, match_ticks_regular, verify_axis, value_interval


def calibrate_figs11(image: np.ndarray):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=130.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    # Note: hl[1] is the horizontal max_rho = 1 homogeneous branch!
    h_one = hl[1] if len(hl) >= 3 else None

    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(6, 22), margin=15, max_width=14)
    # 7 interior x-ticks: 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00 + right border = 2.05
    xt_int = match_ticks_regular(xt_top, 7, tol=3.5)
    x_marks = list(zip(xt_int, [1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00])) + [(v_right[0], 2.05)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.5)

    # 5 interior y-ticks: 10, 8, 6, 4, 2 + homogeneous line max_rho = 1
    yt_int = match_ticks_regular([t for t in yt_left if t < (h_one[0] - 50 if h_one else h_bot[0] - 100)], 5, tol=4.0)
    y_marks = list(zip(yt_int, [10.0, 8.0, 6.0, 4.0, 2.0]))
    if h_one is not None:
        y_marks.append((h_one[0], 1.0))
    y_axis = verify_axis("linear", y_marks, max_residual_px=3.0)
    return x_axis, y_axis, (v_left[0], h_top[0], v_right[0], h_bot[0])


def read(image, at):
    x_axis, y_axis, (cl, rt, cr, rb) = calibrate_figs11(image)
    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180.0) & (rgb[:, :, 1] < 85.0) & (rgb[:, :, 2] < 85.0)

    # Read strictly below the max_rho = 1 line (r_one + 5 .. r_one + 24) where only the bottom halves
    # of the red 'B' labels appear and no bifurcating branches intrude (all branches have max_rho >= 1).
    r_one = float(y_axis.to_pixel(1.0))
    r0 = int(round(r_one)) + 5
    r1 = min(int(rb) - 3, int(round(r_one)) + 26)
    b_cols = clusters(np.where(red[r0:r1, int(cl) + 20 : int(cr) - 20].any(axis=0))[0] + int(cl) + 20, gap=12)
    if len(b_cols) != 7:
        raise RuntimeError(f"Fig. S11: expected 7 bifurcation 'B' markers on max_rho=1, found {len(b_cols)}: {b_cols}")

    # Ascending chi_i^(0) order for alpha=35, phi=0.1 corresponds to modes i = 6, 7, 5, 8, 9, 4, 10
    mode_order = [6, 7, 5, 8, 9, 4, 10]
    out = []
    for mode_i, (c_mid, c_s, c_e) in zip(mode_order, b_cols):
        width = float(c_e - c_s + 1)
        half_px = max(0.5 * width, 4.0) + x_axis.residual_px + 1.0
        val, lo, hi = value_interval(x_axis, c_mid, half_px)
        out.append({
            "x": float(mode_i),
            "y": val,
            "lo": lo,
            "hi": hi,
            "px": [c_mid, r_one],
            "note": f"red 'B' marker on max_rho=1 for mode i={mode_i} (w={width:.0f}px); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })
    return sorted(out, key=lambda r: r["x"])
