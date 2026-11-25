"""Shared calibrated reader for SI Fig. S6 (harmonic spring constant k/F_0 vs L_D)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_box_spines,
    find_ticks,
    match_ticks_regular,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate(image: np.ndarray):
    dark, gray, hs, vs = find_box_spines(image, thresh=140, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # Bottom major ticks at 3, 4, 5; left = 2, right = 6.
    xt_raw = find_ticks(dark, bot, (left[2] + 5, right[1] - 5), "x", side=-1, band=(18, 32), frac=0.5)
    xt = match_ticks_regular(xt_raw, 3, tol=4.0)
    x_marks = [(left[0], 2.0), (xt[0], 3.0), (xt[1], 4.0), (xt[2], 5.0), (right[0], 6.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Right spine has all 7 clean major ticks at 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2 (note bot spine ~ 0.112).
    yt_raw = find_ticks(dark, right, (top[2] + 5, bot[1] - 5), "y", side=-1, band=(18, 32), frac=0.5)
    yt = match_ticks_regular(yt_raw, 7, tol=4.0)
    y_marks = [(yt[i], 1.4 - 0.2 * i) for i in range(7)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    return x_axis, y_axis, (left, top, right, bot)


def read_curve(image: np.ndarray, at: list[float], curve: str) -> list[dict]:
    x_axis, y_axis, (left, top, right, bot) = _calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = image.shape[:2]

    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 3 : bot[1] + 2, left[2] + 3 : right[1] - 3] = True
    # Mask out legend (x in [3.9, 5.5], y in [0.65, 1.15])
    roi[int(y_axis.to_pixel(1.15)) : int(y_axis.to_pixel(0.65)), int(x_axis.to_pixel(3.9)) : int(x_axis.to_pixel(5.6))] = False

    if curve == "LM1p5":
        mask = (B > 180) & (R < 60) & (G < 60) & roi
    elif curve == "LM0p33":
        mask = (R > 180) & (G < 60) & (B < 60) & roi
    else:
        raise ValueError(f"Unknown curve {curve!r}")

    return [sample_solid_curve_y(mask, x_axis, y_axis, float(x), top[2] + 3, bot[1] + 2, search_dx=2) for x in at]
