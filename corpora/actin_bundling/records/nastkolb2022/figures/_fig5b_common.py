"""Shared calibration and dual-axis cross-section curve reader for Fig. 5B (PDF p. 8, 1200 dpi)."""

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


def read_fig5b(image: np.ndarray, channel: str, at: list[float]) -> list[dict]:
    """Read Fig. 5B dual-axis bundle cross-section ('lipids' left blue axis 0.90-1.30, 'actin' right red axis 0.60-1.00)."""
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]
    vr = v_spines[-1]

    # X axis: 0, 2, 4, 6, 8 um (5 major ticks; exclude ticks outside spines)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=42, min_len=14, coord_min=vl[0] - 5, coord_max=vr[0] - 5)
    xt = select_major_ticks(xt_raw, expected=5, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 2.0, 4.0, 6.0, 8.0], max_res_px=2.5)

    if channel == "lipids":
        # Left y-axis: 1.30, 1.25, 1.20, [1.15 tick is missing in PDF!], 1.10, 1.05, 1.00, 0.95, 0.90
        yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=42, min_len=14, coord_min=50, coord_max=hr[1] + 10)
        longest = max(l for _, l in yt_raw)
        yt = [pos for pos, l in yt_raw if l >= 0.72 * longest]
        if len(yt) == 8:
            y_cal = LinearCal(yt, [1.30, 1.25, 1.20, 1.10, 1.05, 1.00, 0.95, 0.90], max_res_px=2.5)
        elif len(yt) == 9:
            y_cal = LinearCal(yt, [1.30, 1.25, 1.20, 1.15, 1.10, 1.05, 1.00, 0.95, 0.90], max_res_px=2.5)
        else:
            raise RuntimeError(f"Unexpected Fig 5B left y ticks: {yt_raw}")
        # Cyan-blue curve (~RGB 64,176,240)
        mask = (B > 175) & (G > 130) & (R < 130) & (B - R > 80)
        # Mask out the "Lipids" text label past x = 8.3 um
        mask[:, int(x_cal.to_px(8.30)) :] = False
    elif channel == "actin":
        # Right y-axis: 1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60 (9 major ticks)
        yt_raw = detect_ticks_along_spine(dark, vr, "y", "right", min_offset=2, max_offset=42, min_len=14, coord_min=50, coord_max=hr[1] + 10)
        yt = select_major_ticks(yt_raw, expected=9, frac=0.72)
        y_cal = LinearCal(yt, [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60], max_res_px=2.5)
        # Red curve (~RGB 245,0,0)
        mask = (R > 180) & (R - G > 110) & (R - B > 110)
        mask[:, int(x_cal.to_px(8.30)) :] = False
    else:
        raise ValueError(channel)

    r_top = 50
    r_bot = hr[0] - 2
    c_left = vl[1] + 1
    c_right = int(x_cal.to_px(8.30))

    return read_colored_curve_at(
        mask, x_cal, y_cal, at, r_top, r_bot, c_left, c_right,
        search_cols=8, window_half=2,
    )
