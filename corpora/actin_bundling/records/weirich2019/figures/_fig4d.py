"""Shared calibrated reader for Main Fig. 4D (spherical droplet division phase diagram, Eq. 7)."""

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

    # 24 evenly spaced x-ticks at 1.2, 1.4, ..., 5.8; left = 1.0, right = 6.0.
    xt_raw = find_ticks(dark, bot, (left[2] + 5, right[1] - 5), "x", side=-1, band=(12, 26), frac=0.5)
    xt = match_ticks_regular(xt_raw, 24, tol=4.0)
    x_marks = [(left[0], 1.0)] + [(xt[i], 1.2 + 0.2 * i) for i in range(24)] + [(right[0], 6.0)]
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

    # Interior region inside spines
    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 4 : bot[1] - 4, left[2] + 4 : right[1] - 4] = True

    if boundary == "unequal":
        # Thick red boundary dots (radius ~14 px); thin red cartoon outlines (< 5 px) killed by 5 erosions.
        red = (R > 160) & ((R - G) > 80) & ((R - B) > 80) & roi
        core = ndimage.binary_erosion(red, iterations=5)
        dots = ndimage.binary_dilation(core, iterations=6) & red
        out = []
        for x in at:
            out.append(
                sample_dashed_curve_y(
                    dots, x_axis, y_axis, float(x), top[2] + 4, bot[1] - 4,
                    max_search_dx=95, sigma_dx=40.0, deg=2, note_prefix="red dots, ",
                )
            )
        return out

    elif boundary == "equal":
        # Black boundary dots (radius ~14 px); exclude narrow axis ticks (width ~8 px) and 3 huge cartoon colloid disks (R~60 px).
        blk = (R < 50) & (G < 50) & (B < 50) & roi
        big_core = ndimage.binary_erosion(blk, iterations=25)
        big = ndimage.binary_dilation(big_core, iterations=40)
        core = ndimage.binary_erosion(blk & (~big), iterations=5)
        dots = ndimage.binary_dilation(core, iterations=6) & blk & (~big)
        out = []
        for x in at:
            out.append(
                sample_dashed_curve_y(
                    dots, x_axis, y_axis, float(x), top[2] + 4, bot[1] - 4,
                    max_search_dx=95, sigma_dx=40.0, deg=2, note_prefix="black dots, ",
                )
            )
        return out
    else:
        raise ValueError(f"Unknown boundary {boundary!r}")
