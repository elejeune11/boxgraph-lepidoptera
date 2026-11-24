"""Shared calibration and curve reader for Fig. 2C transverse bundle intensity profiles (PDF p. 4, 1200 dpi)."""

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


def read_fig2c(image: np.ndarray, time_label: str, at: list[float]) -> list[dict]:
    """Read cleanly separated endpoint profiles in Fig. 2C ('0s' darkest purple, '405s' brightest yellow)."""
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # Legend box right edge / left edge (vertical lines past c=2200)
    leg_c0 = min(c[0] for c in v_spines if c[0] > vl[1] + 1500)

    # X axis: 0, 1, 2, 3, 4 um (5 unobstructed major ticks before the legend border at ~x=4.9)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=45, min_len=14, coord_min=vl[1] + 50, coord_max=leg_c0 - 30)
    xt = select_major_ticks(xt_raw, expected=5, frac=0.70)
    x_cal = LinearCal(xt, [0.0, 1.0, 2.0, 3.0, 4.0], max_res_px=3.0)

    # Y axis: 5e6, 4e6, 3e6, 2e6, 1e6, 0 (6 major ticks)
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=45, min_len=10, coord_min=20, coord_max=hr[1] + 10)
    yt = select_major_ticks(yt_raw, expected=6, frac=0.65)
    y_cal = LinearCal(yt, [5e6, 4e6, 3e6, 2e6, 1e6, 0.0], max_res_px=3.0)

    if time_label == "0s":
        # Darkest purple (bottom curve): R~125, G~0..30, B~125 (R and B similar, G very low)
        mask = (R > 85) & (R < 170) & (B > 85) & (B < 170) & (G < 38) & (np.abs(R - B) < 40)
    elif time_label == "405s":
        # Brightest yellow (top layer): R > 225, G > 220, B < 55
        mask = (R > 220) & (G > 215) & (B < 65)
    else:
        raise ValueError(time_label)

    # Exclude legend interior
    mask[:, leg_c0 - 5 :] = False

    r_top = int(y_cal.to_px(5.05e6))
    r_bot = hr[0] - 4
    c_left = vl[1] + 4
    c_right = leg_c0 - 6

    return read_colored_curve_at(
        mask, x_cal, y_cal, at, r_top, r_bot, c_left, c_right,
        search_cols=12, window_half=3,
    )
