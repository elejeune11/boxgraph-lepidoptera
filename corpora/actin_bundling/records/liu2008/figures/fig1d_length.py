"""Calibrated reader for Fig. 1d (p. 2 / printed p. 790): thin actin protrusion length vs time."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _common as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    box = C.find_plot_box(image, dark_thr=120.0, sat_thr=35.0, min_h_frac=0.55, min_v_frac=0.55)
    dark = C.dark_neutral_mask(image, dark_thr=120.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    # Bottom inward ticks: 0, 5, 10, 15, 20, 25 min (note: curve approaches bottom spine near t = 0)
    xt_all = C.scan_ticks(dark, "x", bot[0] - 26, bot[0] - 6, left[1] + 5, right[0] - 5, min_frac=0.6, max_len=60)
    yt_all = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, left[1] + 6, left[1] + 26, min_frac=0.6, max_len=25)

    # Use the 5 unobstructed x ticks (5, 10, 15, 20, 25) plus verify t=0 tick
    xt_clean = C.scan_ticks(dark, "x", bot[0] - 26, bot[0] - 6, left[1] + 5, right[0] - 5, min_frac=0.6, max_len=15)
    xt_reg = C.select_regular_ticks(xt_clean, 5, tol=2.0)
    xa = C.make_axis("linear", list(zip(xt_reg, [5.0, 10.0, 15.0, 20.0, 25.0])), max_res_px=1.5)

    yt_reg = C.select_regular_ticks(yt_all, 9, tol=2.0)
    ya = C.make_axis("linear", list(zip(yt_reg, [24.0, 22.0, 20.0, 18.0, 16.0, 14.0, 12.0, 10.0, 8.0])), max_res_px=1.5)

    # Thick black elongation polyline inside box, excluding inset kymograph (lower-right: t > 11.5 min, L < 18.5 um)
    black = C.pure_black_mask(image, max_ch=70.0)
    mask = np.zeros_like(black)
    mask[top[1] + 8 : bot[0] - 8, left[1] + 8 : right[0] - 8] = black[top[1] + 8 : bot[0] - 8, left[1] + 8 : right[0] - 8]
    # Remove inward ticks along left and bottom spines (staying clear of curve near t=0, L~7.5)
    mask[:, : left[1] + 30] = False
    mask[bot[0] - 28 :, int(round(xa.to_pixel(2.0))) :] = False
    # Remove inset kymograph in lower right
    mask[int(round(ya.to_pixel(18.5))) :, int(round(xa.to_pixel(11.5))) :] = False

    return C.trace_curve_at(
        mask=mask,
        x_axis=xa,
        y_axis=ya,
        at=at,
        r_min=top[1] + 8,
        r_max=bot[0] - 8,
        half_col=2,
        search_window_px=6,
        min_run=10,
        clamp_min=0.0,
        label="Fig. 1d Length(t)",
    )
