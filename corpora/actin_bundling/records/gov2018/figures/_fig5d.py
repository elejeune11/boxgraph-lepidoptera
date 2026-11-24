"""Shared calibrated reader for Fig. 5(d) (p. 9): steady-state two-curvature wave h(x), n_-(x), n_+(x)."""

from __future__ import annotations

import numpy as np
from _calib import (
    find_inward_ticks,
    find_neutral_frame,
    match_ticks_regular,
    sample_solid_curve_y,
    verify_axis,
)


def _calibrate_fig5d(image: np.ndarray):
    dark, gray, h_lines, v_lines = find_neutral_frame(
        image, row_frac=0.55, col_frac=0.55, gray_thresh=170.0, sat_thresh=35
    )
    if len(v_lines) < 2 or len(h_lines) < 4:
        raise RuntimeError(f"Fig. 5(d): unexpected frame lines h={h_lines}, v={v_lines}")

    v_left, v_right = v_lines[0], v_lines[-1]
    h0, h1, h2, h3 = h_lines[0], h_lines[1], h_lines[2], h_lines[3]
    c0, c1 = int(round(v_left[0])), int(round(v_right[0]))

    # --- Shared x-axis: x = -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10 um ---
    _, xt_b, _, _ = find_inward_ticks(
        dark, int(round(h2[0])), int(round(h3[0])), c0, c1, band=(8, 24), margin=15, max_width=18, frac=0.50
    )
    x_ticks = match_ticks_regular([v_left[0]] + xt_b + [v_right[0]], 11, tol=8.0)
    x_ax = verify_axis("linear", list(zip(x_ticks, [-10.0 + 2.0 * k for k in range(11)])), max_residual_px=5.0)

    # --- Top panel h(x): 1.0 (top border), 0.5, 0.0, -0.5 (left ticks), -1.0 (bottom border) ---
    _, _, yt_h_l, _ = find_inward_ticks(
        dark, int(round(h0[0])), int(round(h1[0])), c0, c1, band=(8, 24), margin=15, max_width=18, frac=0.50
    )
    h_ticks = match_ticks_regular([h0[0]] + yt_h_l + [h1[0]], 5, tol=8.0)
    y_h = verify_axis("linear", list(zip(h_ticks, [1.0, 0.5, 0.0, -0.5, -1.0])), max_residual_px=5.0)

    # --- Middle panel n_-(x): 0.815 (top border), 0.810, 0.805, 0.800, 0.795 (right ticks), 0.790 (bottom border) ---
    _, _, _, yt_m_r = find_inward_ticks(
        dark, int(round(h1[0])), int(round(h2[0])), c0, c1, band=(8, 24), margin=15, max_width=18, frac=0.50
    )
    m_ticks = match_ticks_regular([h1[0]] + yt_m_r + [h2[0]], 6, tol=8.0)
    y_nm = verify_axis("linear", list(zip(m_ticks, [0.815, 0.810, 0.805, 0.800, 0.795, 0.790])), max_residual_px=5.0)

    # --- Bottom panel n_+(x): 0.16 (top border), 0.15, 0.14, 0.13, 0.12 (right ticks), 0.11 (bottom border) ---
    _, _, _, yt_p_r = find_inward_ticks(
        dark, int(round(h2[0])), int(round(h3[0])), c0, c1, band=(8, 24), margin=15, max_width=18, frac=0.50
    )
    p_ticks = match_ticks_regular([h2[0]] + yt_p_r + [h3[0]], 6, tol=8.0)
    y_np = verify_axis("linear", list(zip(p_ticks, [0.16, 0.15, 0.14, 0.13, 0.12, 0.11])), max_residual_px=5.0)

    boxes = {
        "h": (int(h0[2]) + 3, int(h1[1]) - 3),
        "n_convex": (int(h1[2]) + 3, int(h2[1]) - 3),
        "n_concave": (int(h2[2]) + 3, int(h3[1]) - 3),
    }
    return x_ax, {"h": y_h, "n_convex": y_nm, "n_concave": y_np}, boxes


def read_fig5d(image: np.ndarray, at: list[float], panel: str) -> list[dict]:
    x_ax, y_axes, boxes = _calibrate_fig5d(image)
    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    # Dark-blue wave curve (~RGB [101, 116, 179]); excludes black arrow/text and orange concave icon
    mask = (b > 125) & (b - r > 40) & (b - g > 22)
    rmin, rmax = boxes[panel]
    y_ax = y_axes[panel]
    return [
        sample_solid_curve_y(
            mask, x_ax, y_ax, x, rmin, rmax,
            search_dx=3,
            note_prefix=f"fig5d {panel}; ",
        )
        for x in at
    ]
