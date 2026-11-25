"""Shared calibrated reader for SI Fig. S9 (nearly spherical two-tactoid division free energy, Eq. 7)."""

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

    # 21 ticks at 0.00, 0.05, ..., 1.00 (padded inside spines; rasterized bitmap has ~2.5 px quantization).
    xt_raw = find_ticks(dark, bot, (left[2] + 5, right[1] - 5), "x", side=-1, band=(12, 26), frac=0.5)
    xt = match_ticks_regular(xt_raw, 21, tol=8.0)
    x_marks = [(xt[i], 0.05 * i) for i in range(21)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=4.0)

    # 4 major y-ticks on left spine at 2.4, 2.3, 2.2, 2.1.
    yt_raw = find_ticks(dark, left, (top[2] + 5, bot[1] - 5), "y", side=+1, band=(18, 32), frac=0.5)
    yt = match_ticks_regular(yt_raw, 4, tol=8.0)
    y_marks = [(yt[i], 2.4 - 0.1 * i) for i in range(4)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=4.0)

    return x_axis, y_axis, (left, top, right, bot)


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    x_axis, y_axis, (left, top, right, bot) = _calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = image.shape[:2]

    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 35 : bot[1] - 35, left[2] + 45 : right[1] - 45] = True
    # Mask out center legend (x in [0.30, 0.72], y in [2.145, 2.305])
    roi[int(y_axis.to_pixel(2.305)) : int(y_axis.to_pixel(2.145)), int(x_axis.to_pixel(0.30)) : int(x_axis.to_pixel(0.72))] = False

    if curve == "w1p5":
        # Solid purple curve (R ~ 120, G ~ 0, B ~ 255)
        mask = (B > 170) & (R > 60) & (G < 80) & roi
        return [sample_solid_curve_y(mask, x_axis, y_axis, float(x), top[2] + 35, bot[1] - 35, search_dx=3) for x in at]
    elif curve == "w2p5":
        # Dashed red curve
        mask = (R > 180) & (G < 80) & (B < 80) & roi
        return [
            sample_dashed_curve_y(
                mask, x_axis, y_axis, float(x), top[2] + 35, bot[1] - 35,
                max_search_dx=140, sigma_dx=55.0, deg=2, note_prefix="red dashed, ",
            )
            for x in at
        ]
    elif curve == "w3p5":
        # Solid black curve
        mask = (R < 50) & (G < 50) & (B < 50) & roi
        return [sample_solid_curve_y(mask, x_axis, y_axis, float(x), top[2] + 35, bot[1] - 35, search_dx=3) for x in at]
    else:
        raise ValueError(f"Unknown curve {curve!r}")
