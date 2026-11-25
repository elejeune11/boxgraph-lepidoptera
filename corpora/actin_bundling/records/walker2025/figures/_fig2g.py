"""Shared calibrated reader for Fig. 2G: mini-Lpd and VASP FRAP Recovery (%) vs Time (mean +/- SD, n=9)."""

from __future__ import annotations

import numpy as np
import _common as C
from fig1f_frap_minilpd import _read_frap_sd_at_col


def read_fig2g(image: np.ndarray, channel: str, at: list[float]) -> list[dict]:
    ax = C.find_l_axes(image, dark_thresh=110.0, min_h_frac=0.5, min_v_frac=0.5, max_tick_len=40)
    xt = C.match_regular_ticks(ax["x_ticks"], 5, tol=4.0)
    yt = C.match_regular_ticks(ax["y_ticks"], 4, tol=4.0)
    xa = C.make_axis("linear", list(zip(xt, [-10.0, 40.0, 90.0, 140.0, 190.0])), max_res_px=3.5)
    ya = C.make_axis("linear", list(zip(yt, [120.0, 80.0, 40.0, 0.0])), max_res_px=3.5)

    rgb = image.astype(np.float64)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)

    if channel == "minilpd":
        # Emerald green curve: G strongly exceeds R and B
        curve = (G - R > 45.0) & (G - B > 15.0) & (G > 90.0)
    elif channel == "vasp":
        # Bright magenta curve: R and B both high, G low
        curve = (R - G > 70.0) & (B - G > 70.0) & (R > 140.0) & (B > 140.0)
    else:
        raise ValueError(f"Unknown channel: {channel}")

    black = (gray < 100.0) & (sat < 25.0)
    # Exclude upper-right legend area
    c_leg = int(round(xa.to_pixel(115.0)))
    r_leg = int(round(ya.to_pixel(88.0)))
    curve[:r_leg, c_leg:] = False

    r_top = int(round(ya.to_pixel(115.0)))
    r_bot = ax["bot"]["lo"] - 2

    out = []
    for xv in at:
        pc = xa.to_pixel(float(xv))
        out.append(
            _read_frap_sd_at_col(
                image,
                pc,
                curve,
                black,
                r_top,
                r_bot,
                ya,
                float(xv),
                xa,
                note_prefix=f"Fig. 2G {channel} FRAP; ",
            )
        )
    return out
