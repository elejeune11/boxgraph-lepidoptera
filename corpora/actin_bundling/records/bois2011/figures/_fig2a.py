"""Calibration and curve isolation for Fig. 2(a): potentials ell^2 V(u) for Pe = 10, 7, 3."""

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


def calibrate_fig2a(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(dark, h_top, h_bot, v_left, v_right)
    # x ticks: left border (-8), 3 interior ticks (-4, 0, 4), right border (8)
    x_inner = match_ticks_regular(xt_bot, 3, tol=2.0)
    x_axis = verify_axis(
        "linear",
        [(v_left[0], -8.0), (x_inner[0], -4.0), (x_inner[1], 0.0), (x_inner[2], 4.0), (v_right[0], 8.0)],
        max_residual_px=1.0,
    )
    # y ticks: 5 interior ticks on right border (10, 5, 0, -5, -10)
    y_inner = match_ticks_regular(yt_right, 5, tol=2.0)
    y_axis = verify_axis(
        "linear",
        list(zip(y_inner, [10.0, 5.0, 0.0, -5.0, -10.0])),
        max_residual_px=1.0,
    )
    box = (v_left[2] + 15, h_top[2] + 15, v_right[1] - 15, h_bot[1] - 15)
    return dark, x_axis, y_axis, box


def clean_interior_mask(dark: np.ndarray, x_axis, y_axis, box):
    """Mask out the plot borders, top-left legend box, and bottom-center inset box."""
    c0, r0, c1, r1 = box
    mask = np.zeros_like(dark)
    mask[r0 : r1 + 1, c0 : c1 + 1] = dark[r0 : r1 + 1, c0 : c1 + 1]

    # Detect legend box rectangle in upper-left quadrant (u in [-7.8, -0.8], V in [3.0, 12.0])
    lc0 = max(c0, int(x_axis.to_pixel(-7.8)))
    lc1 = min(c1, int(x_axis.to_pixel(-0.8)))
    lr0 = max(r0, int(y_axis.to_pixel(12.0)))
    lr1 = min(r1, int(y_axis.to_pixel(3.0)))
    sub_l = mask[lr0:lr1, lc0:lc1]
    lh = np.where(sub_l.mean(axis=1) > 0.45)[0] + lr0
    lv = np.where(sub_l.mean(axis=0) > 0.45)[0] + lc0
    if len(lh) and len(lv):
        mask[lh.min() - 6 : lh.max() + 6, lv.min() - 6 : lv.max() + 6] = False

    # Detect inset box rectangle in lower-center region (u in [-3.3, 3.5], V in [-13.0, -1.8])
    ic0 = max(c0, int(x_axis.to_pixel(-3.3)))
    ic1 = min(c1, int(x_axis.to_pixel(3.5)))
    ir0 = max(r0, int(y_axis.to_pixel(-1.8)))
    ir1 = min(r1, int(y_axis.to_pixel(-13.0)))
    sub_i = mask[ir0:ir1, ic0:ic1]
    ih = np.where(sub_i.mean(axis=1) > 0.35)[0] + ir0
    iv = np.where(sub_i.mean(axis=0) > 0.35)[0] + ic0
    if len(ih) and len(iv):
        # Also clear inset tick labels (to left and below inset frame)
        mask[ih.min() - 60 : r1 + 1, iv.min() - 110 : iv.max() + 50] = False

    return mask


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read one of 'Pe10' (solid), 'Pe7' (dashed), 'Pe3' (dotted) from Fig. 2(a)."""
    dark, x_axis, y_axis, box = calibrate_fig2a(image)
    mask = clean_interior_mask(dark, x_axis, y_axis, box)
    c0, r0, c1, r1 = box

    # Classify connected components by area:
    # - Pe=10 solid curve is one large component (area > 5000 px)
    # - Pe=7 dashed curve has dashes of area ~200..900 px (and width > height or elongated)
    # - Pe=3 dotted curve has small vertical-ish ticks/dots of area ~30..180 px
    lab, n = ndimage.label(mask)
    counts = np.bincount(lab.ravel())
    solid_mask = np.zeros_like(mask)
    dash_mask = np.zeros_like(mask)
    dot_mask = np.zeros_like(mask)
    for idx in range(1, n + 1):
        area = int(counts[idx])
        if area >= 3000:
            solid_mask |= lab == idx
        elif area >= 190:
            dash_mask |= lab == idx
        elif area >= 25:
            dot_mask |= lab == idx

    # For solid Pe=10: note near u=0 it may briefly touch a dash/dot, which stays in solid_mask;
    # at any x in `at` away from u=0, solid_mask has a single clean run (topmost if touching anything).
    out = []
    for x in at:
        if curve == "Pe10":
            res = sample_solid_curve_y(
                solid_mask,
                x_axis,
                y_axis,
                float(x),
                r0,
                r1,
                pick_fn=lambda rs, col: rs[0],  # Pe=10 is strictly the topmost curve outside [-1.8, 0]
                note_prefix="Pe=10 solid, ",
            )
        elif curve == "Pe7":
            res = sample_dashed_curve_y(
                dash_mask, x_axis, y_axis, float(x), r0, r1, max_search_dx=50, note_prefix="Pe=7 dashed, "
            )
        elif curve == "Pe3":
            res = sample_dashed_curve_y(
                dot_mask, x_axis, y_axis, float(x), r0, r1, max_search_dx=35, note_prefix="Pe=3 dotted, "
            )
        else:
            raise ValueError(f"Unknown Fig. 2(a) curve {curve!r}")
        out.append(res)
    return out
