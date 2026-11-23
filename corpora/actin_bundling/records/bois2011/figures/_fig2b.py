"""Calibration and curve isolation for Fig. 2(b): Pe = 25 potentials for 2-peak and 1-peak steady states."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from _calib import (
    find_frame,
    find_inward_ticks,
    match_ticks_regular,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)


def calibrate_fig2b(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(dark, h_top, h_bot, v_left, v_right)
    # x ticks: left border (-10), 2 interior ticks (-5, 0), right border (5)
    x_inner = match_ticks_regular(xt_bot, 2, tol=2.0)
    x_axis = verify_axis(
        "linear",
        [(v_left[0], -10.0), (x_inner[0], -5.0), (x_inner[1], 0.0), (v_right[0], 5.0)],
        max_residual_px=1.0,
    )
    # y ticks: 5 interior ticks on right border (40, 30, 20, 10, 0)
    y_inner = match_ticks_regular(yt_right, 5, tol=2.0)
    y_axis = verify_axis(
        "linear",
        list(zip(y_inner, [40.0, 30.0, 20.0, 10.0, 0.0])),
        max_residual_px=1.0,
    )
    box = (v_left[2] + 15, h_top[2] + 15, v_right[1] - 15, h_bot[1] - 15)
    return dark, x_axis, y_axis, box


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read '2peak' (solid) or '1peak' (dashed) curve from Fig. 2(b)."""
    dark, x_axis, y_axis, box = calibrate_fig2b(image)
    c0, r0, c1, r1 = box
    mask = np.zeros_like(dark)
    mask[r0 : r1 + 1, c0 : c1 + 1] = dark[r0 : r1 + 1, c0 : c1 + 1]

    # Detect and clear top-center legend box (u in [-7.0, 2.0], V in [30.0, 45.0])
    lc0 = max(c0, int(x_axis.to_pixel(-7.0)))
    lc1 = min(c1, int(x_axis.to_pixel(2.0)))
    lr0 = max(r0, int(y_axis.to_pixel(45.0)))
    lr1 = min(r1, int(y_axis.to_pixel(30.0)))
    sub_l = mask[lr0:lr1, lc0:lc1]
    lh = np.where(sub_l.mean(axis=1) > 0.45)[0] + lr0
    lv = np.where(sub_l.mean(axis=0) > 0.45)[0] + lc0
    if len(lh) and len(lv):
        mask[lh.min() - 6 : lh.max() + 6, lv.min() - 6 : lv.max() + 6] = False

    # Separate solid (2-peak) from dashed (1-peak) by connected component area
    lab, n = ndimage.label(mask)
    counts = np.bincount(lab.ravel())
    solid_mask = np.zeros_like(mask)
    dash_mask = np.zeros_like(mask)
    for idx in range(1, n + 1):
        area = int(counts[idx])
        if area >= 3000:
            solid_mask |= lab == idx
        elif area >= 80:
            dash_mask |= lab == idx

    out = []
    for x in at:
        if curve == "2peak":
            # At u < 0, solid is upper curve; at u > 0, solid is lower curve (if a dash touches near u=0)
            picker = (lambda rs, col: rs[0]) if float(x) < 0 else (lambda rs, col: rs[-1])
            res = sample_solid_curve_y(
                solid_mask, x_axis, y_axis, float(x), r0, r1, pick_fn=picker, note_prefix="2-peak solid, "
            )
        elif curve == "1peak":
            res = sample_dashed_curve_y(
                dash_mask, x_axis, y_axis, float(x), r0, r1, max_search_dx=60, note_prefix="1-peak dashed, "
            )
        else:
            raise ValueError(f"Unknown Fig. 2(b) curve {curve!r}")
        out.append(res)
    return out
