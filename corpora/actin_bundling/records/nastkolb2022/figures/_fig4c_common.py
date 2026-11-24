"""Shared calibration and curve reader for Fig. 4C FRAP recovery curves (PDF p. 7, 1200 dpi)."""

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


def read_fig4c(image: np.ndarray, curve: str, at: list[float]) -> list[dict]:
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: -0.2 (at y-axis spine corner), 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 min
    # Fit on the 7 unobstructed ticks 0.0..1.2 (exact 444.0 px spacing, 0.00 px residual)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=55, min_len=14, coord_min=vl[0] - 10, coord_max=W - 20)
    xt = select_major_ticks(xt_raw, expected=8, frac=0.72)
    x_cal = LinearCal(xt[1:], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2], max_res_px=1.5)

    # Y axis: 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=46, min_len=14, coord_min=50, coord_max=hr[0] - 5)
    yt = select_major_ticks(yt_raw, expected=7, frac=0.72)
    y_cal = LinearCal(yt, [1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0], max_res_px=2.0)

    # Mask out legend box (x in [0.70, 1.38] min, y in [-0.05, 0.42]) so legend strokes never interfere!
    leg_c0 = int(x_cal.to_px(0.70))
    leg_c1 = W
    leg_r0 = int(y_cal.to_px(0.42))
    leg_r1 = hr[0]

    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    if curve == "slb_prepol":
        # Dark gray curve (~RGB 81,81,81)
        mask = (gray > 45) & (gray < 125) & (sat < 22)
    elif curve == "slb_pol30":
        # Red curve (~RGB 241,64,64)
        mask = (R > 180) & (R - G > 110) & (R - B > 110)
    elif curve == "ml_prepol":
        # Blue curve (~RGB 26,111,223)
        mask = (B > 155) & (B - R > 110) & (B - G > 45)
    elif curve == "ml_pol30":
        # Green curve (~RGB 55,173,107)
        mask = (G > 125) & (G - R > 60) & (G - B > 25)
    else:
        raise ValueError(curve)

    mask[leg_r0:leg_r1, leg_c0:leg_c1] = False

    r_top = int(y_cal.to_px(1.15))
    r_bot = hr[0] - 4
    c_left = vl[1] + 4
    c_right = W - 10

    # At t = 0.00 min, the vertical photobleach drop line spans y=1.0 down to y=0.0;
    # restrict y_hint at t <= 0.01 to [-0.08, 0.08] (the first post-bleach minimum!)
    def y_hint(xv: float) -> tuple[float, float]:
        if abs(xv) <= 0.01:
            return (-0.08, 0.08)
        return (-0.08, 1.12)

    return read_colored_curve_at(
        mask, x_cal, y_cal, at, r_top, r_bot, c_left, c_right,
        search_cols=10, window_half=3, y_hint_fn=y_hint,
    )
