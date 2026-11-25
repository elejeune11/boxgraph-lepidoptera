"""Shared calibrated reader for Main Fig. 2F (continuum-model centering free energy vs d/L_D)."""

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
    dark, gray, hs, vs = find_box_spines(image, thresh=100, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # Top spine has 3 major ticks at -0.5, 0.0, 0.5 (minor ticks are shorter); left = -1.0, right = +1.0.
    xt_raw = find_ticks(dark, top, (left[2] + 10, right[1] - 10), "x", side=+1, band=(12, 26), frac=0.5)
    xt = match_ticks_regular(xt_raw, 3, tol=5.0)
    x_marks = [(left[0], -1.0), (xt[0], -0.5), (xt[1], 0.0), (xt[2], 0.5), (right[0], 1.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # Left spine has 5 major ticks at 10, 8, 6, 4, 2; bottom spine = 0.0.
    yt_raw = find_ticks(dark, left, (top[2] + 10, bot[1] - 10), "y", side=+1, band=(12, 26), frac=0.5)
    yt = match_ticks_regular(yt_raw, 5, tol=5.0)
    y_marks = [(yt[i], 10.0 - 2.0 * i) for i in range(5)] + [(bot[0], 0.0)]
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
    # Exclude legend box (top-right) and inset tactoid schematic (middle-upper)
    roi[: top[2] + 260, int(x_axis.to_pixel(-0.08)) : int(x_axis.to_pixel(0.62))] = False
    # Exclude tactoid cartoon region (red director lines lie in d/L_D in [-0.45, +0.45], y in [3.2, 7.5])
    r_cartoon_top = int(y_axis.to_pixel(7.6))
    r_cartoon_bot = int(y_axis.to_pixel(3.0))
    c_cartoon_lo = int(x_axis.to_pixel(-0.46))
    c_cartoon_hi = int(x_axis.to_pixel(0.46))

    if curve == "LM1p5":
        mask = (B > 110) & ((B - R) > 40) & ((B - G) > 35) & roi
        row_min, row_max = top[2] + 3, bot[1] + 2
    elif curve == "LM0p33":
        mask = (R > 160) & ((R - G) > 80) & ((R - B) > 80) & roi
        mask[r_cartoon_top:r_cartoon_bot, c_cartoon_lo:c_cartoon_hi] = False
        row_min, row_max = top[2] + 3, bot[1] + 2
    else:
        raise ValueError(f"Unknown curve {curve!r}")

    out = []
    for x in at:
        out.append(sample_solid_curve_y(mask, x_axis, y_axis, float(x), row_min, row_max, search_dx=2))
    return out
