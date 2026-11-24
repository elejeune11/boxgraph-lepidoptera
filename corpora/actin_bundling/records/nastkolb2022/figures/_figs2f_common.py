"""Shared calibration and 5-curve reader for Fig. S2F raw TIRF actin intensity kinetics (PDF p. 16, 1200 dpi)."""

from __future__ import annotations

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (
    LinearCal,
    detect_ticks_along_spine,
    find_black_spines,
    read_colored_curve_at,
    select_major_ticks,
)


def read_figs2f(image: np.ndarray, series: str, at: list[float]) -> list[dict]:
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: 0, 10, 20, 30, 40, 50, 60 min (7 major ticks)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=48, min_len=12, coord_min=vl[1], coord_max=W)
    xt = select_major_ticks(xt_raw, expected=7, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], max_res_px=2.5)

    # Y axis: 10000, 8000, 6000, 4000, 2000, 0 AU (6 major ticks)
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=48, min_len=12, coord_min=50, coord_max=hr[0] - 5)
    yt = select_major_ticks(yt_raw, expected=6, frac=0.72)
    y_cal = LinearCal(yt, [10000.0, 8000.0, 6000.0, 4000.0, 2000.0, 0.0], max_res_px=2.5)

    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    if series == "005":
        mask = (gray > 50) & (gray < 125) & (sat < 22)
    elif series == "0125":
        mask = (R > 180) & (R - G > 110) & (R - B > 110)
    elif series == "025":
        mask = (B > 155) & (B - R > 110) & (B - G > 40)
    elif series == "050":
        mask = (G > 125) & (G - R > 60) & (G - B > 25)
    elif series == "100":
        mask = (R > 140) & (B > 175) & (G < 155) & (B - G > 50) & (R - G > 25)
    else:
        raise ValueError(series)

    # Mask out legend box (x in [34, 60], y in [2200, 5300])
    mask[int(y_cal.to_px(5300.0)) : int(y_cal.to_px(2200.0)), int(x_cal.to_px(34.0)) :] = False

    r_top = int(y_cal.to_px(10100.0))
    r_bot = hr[0] - 3
    c_left = vl[1] + 4
    c_right = W - 2

    return read_colored_curve_at(
        mask, x_cal, y_cal, at, r_top, r_bot, c_left, c_right,
        search_cols=8, window_half=2,
    )
