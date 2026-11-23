"""Calibrated reader for SI Fig. S13 inset (p. 22): contractility Peclet number ramp chi(t) vs time (s)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import clusters, find_inward_ticks, match_ticks_regular, verify_axis, sample_solid_curve_y


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 110.0
    H, W = dark.shape

    # Locate the inset frame strictly inside the main plot box
    r_s0, r_s1 = int(0.22 * H), int(0.80 * H)
    c_s0, c_s1 = int(0.32 * W), int(0.90 * W)
    ih = clusters(np.where(dark[r_s0:r_s1, c_s0:c_s1].mean(axis=1) > 0.60)[0] + r_s0)
    iv = clusters(np.where(dark[r_s0:r_s1, c_s0:c_s1].mean(axis=0) > 0.50)[0] + c_s0)
    if len(ih) < 2 or len(iv) < 2:
        raise RuntimeError(f"Fig. S13 inset: could not locate inset box: ih={ih}, iv={iv}")
    h_top, h_bot = ih[0], ih[-1]
    v_left, v_right = iv[0], iv[-1]

    _, xt_bot, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 16), margin=12, max_width=12)
    # X ticks: left=0, 2 interior = 1000, 2000, right=3000
    xt_int = match_ticks_regular(xt_bot, 2, tol=3.5)
    x_axis = verify_axis("linear", [(v_left[0], 0.0), (xt_int[0], 1000.0), (xt_int[1], 2000.0), (v_right[0], 3000.0)], max_residual_px=2.5)
    # Y ticks: top=1.5, 4 interior = 1.4, 1.3, 1.2, 1.1, bottom=1.0
    yt_int = match_ticks_regular(yt_left, 4, tol=3.5)
    y_axis = verify_axis(
        "linear",
        [(h_top[0], 1.5)] + list(zip(yt_int, [1.4, 1.3, 1.2, 1.1])) + [(h_bot[0], 1.0)],
        max_residual_px=2.5,
    )

    rgb = image.astype(float)
    blue = (rgb[:, :, 2] > 180.0) & (rgb[:, :, 0] < 60.0) & (rgb[:, :, 1] < 60.0)
    r_min, r_max = int(h_top[2]) + 4, int(h_bot[1]) - 4
    return [sample_solid_curve_y(blue, x_axis, y_axis, x, r_min, r_max, search_dx=2) for x in at]
