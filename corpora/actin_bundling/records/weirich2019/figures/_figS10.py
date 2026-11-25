"""Shared calibrated reader for SI Fig. S10 (highly elongated two-tactoid division free energy, Eq. 8)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_box_spines,
    find_ticks,
    match_ticks_regular,
    sample_dashed_curve_y,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate(image: np.ndarray):
    dark, gray, hs, vs = find_box_spines(image, thresh=140, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # 21 evenly spaced x-ticks at 0.00, 0.05, ..., 1.00 (padded inside spines).
    xt_raw = find_ticks(dark, bot, (left[2] + 5, right[1] - 5), "x", side=-1, band=(12, 26), frac=0.5)
    xt = match_ticks_regular(xt_raw, 21, tol=4.0)
    x_marks = [(xt[i], 0.05 * i) for i in range(21)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # 17 evenly spaced y-ticks at 2.1, 2.0, ..., 0.5 (major labeled ticks at 2.0, 1.5, 1.0, 0.5).
    yt_raw = find_ticks(dark, left, (top[2] + 5, bot[1] - 5), "y", side=+1, band=(12, 26), frac=0.5)
    yt = match_ticks_regular(yt_raw, 17, tol=4.0)
    y_marks = [(yt[i], 2.1 - 0.1 * i) for i in range(17)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    x_axis, y_axis, (left, top, right, bot) = _calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = image.shape[:2]

    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 35 : bot[1] - 35, left[2] + 45 : right[1] - 45] = True
    # Mask out center legend (x in [0.24, 0.78], y in [1.00, 1.78])
    roi[int(y_axis.to_pixel(1.78)) : int(y_axis.to_pixel(1.00)), int(x_axis.to_pixel(0.24)) : int(x_axis.to_pixel(0.78))] = False

    if curve == "w0p1":
        # Solid blue curve
        mask = (B > 180) & (R < 60) & (G < 60) & roi
        return [sample_solid_curve_y(mask, x_axis, y_axis, float(x), top[2] + 35, bot[1] - 35, search_dx=3) for x in at]
    elif curve == "w0p5":
        # Dashed red curve
        mask = (R > 180) & (G < 60) & (B < 60) & roi
        return [
            sample_dashed_curve_y(
                mask, x_axis, y_axis, float(x), top[2] + 35, bot[1] - 35,
                max_search_dx=140, sigma_dx=55.0, deg=2, note_prefix="red dashed, ",
            )
            for x in at
        ]
    elif curve == "w2p0":
        # Solid black curve
        mask = (R < 50) & (G < 50) & (B < 50) & roi
        return [sample_solid_curve_y(mask, x_axis, y_axis, float(x), top[2] + 35, bot[1] - 35, search_dx=3) for x in at]
    else:
        raise ValueError(f"Unknown curve {curve!r}")
