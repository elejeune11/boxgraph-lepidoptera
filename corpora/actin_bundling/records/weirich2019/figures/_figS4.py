"""Shared calibrated reader for SI Fig. S4 (dimensionless rod centering free energy vs x_c/L_D)."""

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

    # Bottom major ticks at -0.5, 0.0, 0.5; left = -1.0, right = +1.0.
    xt_raw = find_ticks(dark, bot, (left[2] + 5, right[1] - 5), "x", side=-1, band=(18, 32), frac=0.5)
    xt = match_ticks_regular(xt_raw, 3, tol=4.0)
    x_marks = [(left[0], -1.0), (xt[0], -0.5), (xt[1], 0.0), (xt[2], 0.5), (right[0], 1.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Left major ticks at 1.5, 1.0, 0.5; bottom spine = 0.0.
    yt_raw = find_ticks(dark, left, (top[2] + 5, bot[1] - 5), "y", side=+1, band=(18, 32), frac=0.5)
    yt = match_ticks_regular(yt_raw, 3, tol=4.0)
    y_marks = [(yt[0], 1.5), (yt[1], 1.0), (yt[2], 0.5), (bot[0], 0.0)]
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
    # Mask out center legend (x in [-0.30, +0.35], y in [0.50, 1.35])
    roi[int(y_axis.to_pixel(1.38)) : int(y_axis.to_pixel(0.45)), int(x_axis.to_pixel(-0.32)) : int(x_axis.to_pixel(0.38))] = False

    if curve == "ratio0p4":
        mask = (B > 180) & (R < 60) & (G < 60) & roi
    elif curve == "ratio0p2":
        mask = (R > 180) & (G < 60) & (B < 60) & roi
    elif curve == "ratio0p1":
        mask = (G > 180) & (R < 60) & (B < 60) & roi
    else:
        raise ValueError(f"Unknown curve {curve!r}")

    return [sample_solid_curve_y(mask, x_axis, y_axis, float(x), top[2] + 3, bot[1] + 2, search_dx=2) for x in at]
