"""Shared calibration and curve isolation for Fig. 5(d): ASDM steady-state profiles a, s, v/U."""

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


def calibrate_fig5d(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_frame(image, row_frac=0.5, col_frac=0.45, thresh=100.0)
    h_top, h_bot = h_lines[0], h_lines[-1]
    v_left, v_right = v_lines[0], v_lines[-1]

    xt_top, xt_bot, yt_left, yt_right = find_inward_ticks(dark, h_top, h_bot, v_left, v_right)
    # x ticks: left border (0), 3 interior ticks (5, 10, 15), right border (20)
    x_inner = match_ticks_regular(xt_bot, 3, tol=2.0)
    x_axis = verify_axis(
        "linear",
        [(v_left[0], 0.0), (x_inner[0], 5.0), (x_inner[1], 10.0), (x_inner[2], 15.0), (v_right[0], 20.0)],
        max_residual_px=1.0,
    )
    # left y ticks (a, s): 5 interior ticks at 5, 4, 3, 2, 1
    y_left_inner = match_ticks_regular(yt_left, 5, tol=2.0)
    y_left = verify_axis("linear", list(zip(y_left_inner, [5.0, 4.0, 3.0, 2.0, 1.0])), max_residual_px=1.0)
    # right y ticks (v/U): 3 regular interior ticks at +0.2, 0.0, -0.2
    y_right_inner = match_ticks_regular(yt_right, 3, tol=2.0)
    y_right = verify_axis("linear", list(zip(y_right_inner, [0.2, 0.0, -0.2])), max_residual_px=1.0)

    # Keep r0 right below top border (ht[2]+2) because peaks of a reach very close to top frame
    box = (v_left[2] + 15, h_top[2] + 2, v_right[1] - 15, h_bot[1] - 8)
    return gray, x_axis, y_left, y_right, box


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    """Read 'a' (black solid), 's' (gray solid), or 'v' (black dashed) from Fig. 5(d)."""
    gray, x_axis, y_left, y_right, box = calibrate_fig5d(image)
    c0, r0, c1, r1 = box

    if curve == "s":
        # Gray substrate curve has intensity ~149 (115 <= gray <= 185); black curves are ~33.
        # Opening with a 5x5 disk eliminates 1-2 px antialiased gray fringes along black strokes.
        raw_gray = np.zeros(gray.shape, dtype=bool)
        raw_gray[r0 : r1 + 1, c0 : c1 + 1] = (gray[r0 : r1 + 1, c0 : c1 + 1] >= 115.0) & (
            gray[r0 : r1 + 1, c0 : c1 + 1] <= 185.0
        )
        yy, xx = np.ogrid[-2:3, -2:3]
        disk = (xx * xx + yy * yy) <= 5
        core = ndimage.binary_opening(raw_gray, structure=disk)
        lab, n = ndimage.label(core)
        counts = np.bincount(lab.ravel())
        clean_core = np.zeros_like(core)
        for idx in range(1, n + 1):
            if counts[idx] >= 300:
                clean_core |= lab == idx
        # Restore exact gray stroke edges without picking up black-line fringes
        s_mask = raw_gray & ndimage.binary_dilation(clean_core, structure=np.ones((7, 7), bool))
        return [
            sample_solid_curve_y(s_mask, x_axis, y_left, float(x), r0, r1, note_prefix="Fig. 5(d) s gray, ")
            for x in at
        ]

    # Black curves only (gray < 85 excludes gray substrate curve completely)
    blk = np.zeros(gray.shape, dtype=bool)
    blk[r0 : r1 + 1, c0 : c1 + 1] = gray[r0 : r1 + 1, c0 : c1 + 1] < 85.0

    lab, n = ndimage.label(blk)
    counts = np.bincount(lab.ravel())
    solid_mask = np.zeros_like(blk)
    dash_mask = np.zeros_like(blk)
    for idx in range(1, n + 1):
        area = int(counts[idx])
        if area >= 8000:
            solid_mask |= lab == idx
        elif area >= 80:
            dash_mask |= lab == idx

    if curve == "a":
        return [
            sample_solid_curve_y(solid_mask, x_axis, y_left, float(x), r0, r1, note_prefix="Fig. 5(d) a solid, ")
            for x in at
        ]
    elif curve == "v":
        return [
            sample_dashed_curve_y(
                dash_mask, x_axis, y_right, float(x), r0, r1, max_search_dx=55, note_prefix="Fig. 5(d) v/U dashed, "
            )
            for x in at
        ]
    else:
        raise ValueError(f"Unknown Fig. 5(d) curve {curve!r}")
