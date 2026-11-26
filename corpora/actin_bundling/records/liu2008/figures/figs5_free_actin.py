"""Calibrated reader for Fig. s5 (p. 14 / SI p. 9): normalized free G-actin concentration c_0(t)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    box = C.find_plot_box(image, dark_thr=120.0, sat_thr=35.0, min_h_frac=0.55, min_v_frac=0.65)
    dark = C.dark_neutral_mask(image, dark_thr=120.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    # Top inward x ticks: 400, 800, 1200, 1600 (t=0 is occluded at top by red curve at y=1.0, visible at bottom)
    xt_top = C.scan_ticks(dark, "x", top[1] + 10, top[1] + 40, left[1] + 5, right[0] - 5, min_frac=0.5, min_len=6, max_len=25)
    xt_bot = C.scan_ticks(dark, "x", bot[0] - 40, bot[0] - 10, left[1] + 5, right[0] - 5, min_frac=0.5, min_len=6, max_len=25)
    # Combine unique tick candidates across top and bottom spines
    xt_comb = sorted(xt_top + [x for x in xt_bot if all(abs(x - t) > 25.0 for t in xt_top)])
    xt_reg = C.select_regular_ticks(xt_comb, 5, tol=8.0)
    xa = C.make_axis("linear", list(zip(xt_reg, [0.0, 400.0, 800.0, 1200.0, 1600.0])), max_res_px=3.0)

    # Y ticks: left inward has 0.8, 0.6, 0.4, 0.2, 0.0; right inward has 1.0, 0.8, 0.6, 0.4, 0.2
    yt_l = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, left[1] + 10, left[1] + 42, min_frac=0.5, min_len=6, max_len=25)
    yt_r = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, right[0] - 42, right[0] - 10, min_frac=0.5, min_len=6, max_len=25)
    yt_comb = sorted(yt_r + [y for y in yt_l if all(abs(y - t) > 25.0 for t in yt_r)])
    yt_reg = C.select_regular_ticks(yt_comb, 6, tol=8.0)
    ya = C.make_axis("linear", list(zip(yt_reg, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_res_px=3.5)

    rgb = image.astype(np.float32)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    red = np.zeros(image.shape[:2], dtype=bool)
    red[top[1] + 2 : bot[0] - 2, left[1] + 2 : right[0] - 2] = (
        (R[top[1] + 2 : bot[0] - 2, left[1] + 2 : right[0] - 2] > 140.0)
        & ((R[top[1] + 2 : bot[0] - 2, left[1] + 2 : right[0] - 2] - G[top[1] + 2 : bot[0] - 2, left[1] + 2 : right[0] - 2]) > 75.0)
        & ((R[top[1] + 2 : bot[0] - 2, left[1] + 2 : right[0] - 2] - B[top[1] + 2 : bot[0] - 2, left[1] + 2 : right[0] - 2]) > 75.0)
    )

    return C.trace_curve_at(
        mask=red,
        x_axis=xa,
        y_axis=ya,
        at=at,
        r_min=top[1] + 2,
        r_max=bot[0] - 2,
        half_col=2,
        search_window_px=6,
        min_run=5,
        clamp_min=0.0,
        clamp_max=1.0,
        label="Fig. s5 c_0(t)",
    )
