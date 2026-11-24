"""Shared calibration and dual-axis boxplot reader for Fig. S3F actin titration (PDF p. 17, 1200 dpi)."""

from __future__ import annotations

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (
    LinearCal,
    detect_ticks_along_spine,
    find_black_spines,
    read_origin_boxplots,
    select_major_ticks,
)


def read_figs3f(image: np.ndarray, series: str) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    dark = gray < 105.0
    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]
    vr = v_spines[-1]

    # X axis: 0, 1, 2, 3, 4 uM actin (5 major ticks); plotted groups at 0.5, 1.0, 2.0, 4.0 uM
    x_concs = [0.5, 1.0, 2.0, 4.0]
    r_plot_top = 40
    r_plot_bot = hr[0] - 2
    c_plot_left = vl[1] + 5
    c_plot_right = vr[0] - 5

    if series == "width":
        # Left y-axis: 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0 um (7 major ticks)
        yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=48, min_len=14, coord_min=40, coord_max=hr[0] - 5)
        yt = select_major_ticks(yt_raw, expected=7, frac=0.72)
        y_cal = LinearCal(yt, [3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0], max_res_px=2.5)
        return read_origin_boxplots(
            image, x_concs, "red", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=True, min_box_area=120,
        )
    elif series == "pore":
        # Right y-axis: 12, 10, 8, 6, 4, 2, 0 um^2 (7 major ticks)
        yt_raw = detect_ticks_along_spine(dark, vr, "y", "right", min_offset=2, max_offset=48, min_len=14, coord_min=40, coord_max=hr[0] - 5)
        yt = select_major_ticks(yt_raw, expected=7, frac=0.70)
        y_cal = LinearCal(yt, [12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0], max_res_px=2.5)
        return read_origin_boxplots(
            image, x_concs, "blue", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=True, min_box_area=120,
        )
    else:
        raise ValueError(series)
