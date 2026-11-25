"""Shared calibrated reader for SI Fig. S11 (elongated tactoid division phase diagram, Eq. 8)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from _calib import (
    find_box_spines,
    find_ticks,
    match_ticks_regular,
    match_ticks_values,
    sample_dashed_curve_y,
    verify_axis,
)


def _calibrate(image: np.ndarray):
    dark, gray, hs, vs = find_box_spines(image, thresh=140, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # Bottom major ticks at 2, 3, 4; left = 1, right = 5.
    xt_raw = find_ticks(dark, bot, (left[2] + 5, right[1] - 5), "x", side=-1, band=(18, 32), frac=0.5)
    xt = match_ticks_regular(xt_raw, 3, tol=4.0)
    x_marks = [(left[0], 1.0), (xt[0], 2.0), (xt[1], 3.0), (xt[2], 4.0), (right[0], 5.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Log10 y-axis: 4 major ticks at 20, 10, 5, 2 + bottom spine = 1.
    yt_maj = find_ticks(dark, left, (top[2] + 5, bot[1] - 5), "y", side=+1, band=(18, 32), frac=0.5)
    y_marks = match_ticks_values(yt_maj + [bot[0]], [1.0, 2.0, 5.0, 10.0, 20.0], scale="log10", tol=4.0, increasing_pixel=False)
    y_axis = verify_axis("log10", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def read_boundary(image: np.ndarray, at: list[float], boundary: str) -> list[dict]:
    x_axis, y_axis, (left, top, right, bot) = _calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = image.shape[:2]

    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 4 : bot[1] - 4, left[2] + 4 : right[1] - 4] = True

    if boundary == "unequal":
        # Thick red dots (radius ~18 px at 900 dpi); thin red cartoon outlines (< 5 px) killed by 6 erosions.
        red = (R > 180) & ((R - G) > 100) & ((R - B) > 100) & roi
        core = ndimage.binary_erosion(red, iterations=6)
        dots = ndimage.binary_dilation(core, iterations=7) & red
        return [
            sample_dashed_curve_y(
                dots, x_axis, y_axis, float(x), top[2] + 4, bot[1] - 4,
                max_search_dx=110, sigma_dx=45.0, deg=2, note_prefix="red dots, ",
            )
            for x in at
        ]
    elif boundary == "equal":
        # Black boundary dots (radius ~18 px); exclude axis ticks and 3 huge black cartoon colloid disks.
        blk = (R < 50) & (G < 50) & (B < 50) & roi
        big_core = ndimage.binary_erosion(blk, iterations=32)
        big = ndimage.binary_dilation(big_core, iterations=55)
        core = ndimage.binary_erosion(blk & (~big), iterations=6)
        dots = ndimage.binary_dilation(core, iterations=7) & blk & (~big)
        return [
            sample_dashed_curve_y(
                dots, x_axis, y_axis, float(x), top[2] + 4, bot[1] - 4,
                max_search_dx=110, sigma_dx=45.0, deg=2, note_prefix="black dots, ",
            )
            for x in at
        ]
    else:
        raise ValueError(f"Unknown boundary {boundary!r}")
