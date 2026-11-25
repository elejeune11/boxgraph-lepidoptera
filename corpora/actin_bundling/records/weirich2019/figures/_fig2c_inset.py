"""Shared calibrated reader for Main Fig. 2C Inset (single trajectory d vs Time and exponential fit)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage, signal
from _calib import (
    find_box_spines,
    find_ticks,
    make_point_reading,
    match_ticks_regular,
    sample_dashed_curve_y,
    verify_axis,
)


def _calibrate(image: np.ndarray):
    dark, gray, hs, vs = find_box_spines(image, thresh=100, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # 4 interior x-ticks at 200, 400, 600, 800 s; left = 0, right = 1000.
    xt_raw = find_ticks(dark, top, (left[2] + 10, right[1] - 10), "x", side=+1, band=(5, 25), frac=0.5)
    xt = match_ticks_regular(xt_raw, 4, tol=5.0)
    x_marks = [(left[0], 0.0)] + [(xt[i], 200.0 * (i + 1)) for i in range(4)] + [(right[0], 1000.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Inverted y-axis: top = -3, interior ticks = -2, -1, 0, bot = +1 um.
    yt_raw = find_ticks(dark, left, (top[2] + 10, bot[1] - 10), "y", side=+1, band=(5, 25), frac=0.5)
    yt = match_ticks_regular(yt_raw, 3, tol=5.0)
    y_marks = [(top[0], -3.0)] + [(yt[i], -2.0 + float(i)) for i in range(3)] + [(bot[0], 1.0)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def _separate_red(image: np.ndarray, spines):
    left, top, right, bot = spines
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)

    H, W = image.shape[:2]
    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 4 : bot[1] - 4, left[2] + 4 : right[1] - 4] = True
    red = (R > 150) & ((R - G) > 80) & ((R - B) > 80) & roi

    # Circular dots have radius ~14 px at 2400 dpi; thin curve is ~5 px thick.
    # Iterative 3x3 erosion/dilation is ~20x faster than full-image EDT.
    dots_core = ndimage.binary_erosion(red, iterations=6)
    dots_dil = ndimage.binary_dilation(dots_core, iterations=12)
    curve_only = red & (~dots_dil)
    return red, dots_core, curve_only


def read_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, spines = _calibrate(image)
    left, top, right, bot = spines
    _, _, curve_only = _separate_red(image, spines)
    out = []
    for t in at:
        out.append(
            sample_dashed_curve_y(
                curve_only,
                x_axis,
                y_axis,
                float(t),
                top[2] + 4,
                bot[1] - 4,
                max_search_dx=140,
                sigma_dx=55.0,
                deg=2,
                note_prefix="exp fit, ",
            )
        )
    return out


def read_data(image: np.ndarray) -> list[dict]:
    """Extract individual red trajectory dots by matching pursuit on the opened dot mask."""
    x_axis, y_axis, spines = _calibrate(image)
    left, top, right, bot = spines
    red, dots_core, _ = _separate_red(image, spines)
    dots_mask = ndimage.binary_dilation(dots_core, iterations=7) & red

    # Dot radius at 2400 dpi is ~14 px (diameter ~28 px; adjacent 10-s frames are ~19 px apart in x).
    from _calib import detect_template_mp

    rad = 14.0
    rad_i = 14
    yy, xx = np.ogrid[-rad_i : rad_i + 1, -rad_i : rad_i + 1]
    tpl = ((xx * xx + yy * yy) <= rad * rad).astype(float)
    found = detect_template_mp(dots_mask, tpl, min_full_frac=0.60, min_res_frac=0.45, min_dist=13, refine_r=4)
    found.sort(key=lambda p: p[1])
    cal_note = f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
    return [
        make_point_reading(c, r, x_axis, y_axis, half_px_y=3.5, note=f"dot R=14px, {cal_note}")
        for r, c in found
    ]
