"""Shared calibrated reader for Fig. 4(d) inset top panel (p. 8): h(x) at t = 60, 90, 120 s."""

from __future__ import annotations

import numpy as np
from _calib import (
    clusters,
    find_inward_ticks,
    find_neutral_frame,
    match_ticks_regular,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate_inset(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_neutral_frame(
        image, row_frac=0.55, col_frac=0.55, gray_thresh=155.0, sat_thresh=40
    )
    if len(v_lines) < 2 or len(h_lines) < 2:
        raise RuntimeError(f"Fig. 4(d) inset: unexpected frame lines h={h_lines}, v={v_lines}")

    v_left, v_right = v_lines[0], v_lines[-1]
    h_top, h_bot = h_lines[0], h_lines[1]
    c0, c1 = int(round(v_left[0])), int(round(v_right[0]))
    r0, r1 = int(round(h_top[0])), int(round(h_bot[0]))

    # --- Top x-axis: x = 0, 5, 10, 15, 20, 25 um ---
    xt_t, _, _, _ = find_inward_ticks(
        gray < 175.0, r0, r1, c0, c1, band=(8, 26), margin=15, max_width=28, frac=0.35
    )
    x_ticks = match_ticks_regular([v_left[0]] + xt_t + [v_right[0]], 6, tol=10.0)
    x_ax = verify_axis("linear", list(zip(x_ticks, [0.0, 5.0, 10.0, 15.0, 20.0, 25.0])), max_residual_px=5.0)

    # --- y-axis: 0.05, 0, -0.05, -0.10 um ---
    # Right border has clean ticks at all 4 values (and horizontal zero line confirms 0)
    _, _, _, yt_r = find_inward_ticks(
        gray < 175.0, r0, r1, c0, c1, band=(8, 26), margin=15, max_width=28, frac=0.35
    )
    y_ticks = match_ticks_regular(yt_r, 4, tol=10.0)
    y_ax = verify_axis("linear", list(zip(y_ticks, [0.05, 0.0, -0.05, -0.10])), max_residual_px=5.0)

    box = (int(h_top[2]) + 4, int(h_bot[1]) - 4, int(v_left[2]) + 4, int(v_right[1]) - 4)
    return x_ax, y_ax, box


def _color_masks(image: np.ndarray) -> dict[str, np.ndarray]:
    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    return {
        "t60": (b > 145) & (b - r > 38) & (b - g > 18),
        "t90": (g > 110) & (g - r > 10) & (g - b > 3) & (r > 85) & (b < 175),
        "t120": (r > 170) & (r - g > 42) & (r - b > 42),
    }


def read_fig4d_inset(image: np.ndarray, at: list[float], time_key: str) -> list[dict]:
    x_ax, y_ax, box = _calibrate_inset(image)
    mask = _color_masks(image)[time_key]
    rmin, rmax, _, _ = box
    return [
        sample_solid_curve_y(
            mask, x_ax, y_ax, x, rmin, rmax,
            search_dx=8,
            note_prefix=f"fig4d inset h({time_key}); ",
        )
        for x in at
    ]
