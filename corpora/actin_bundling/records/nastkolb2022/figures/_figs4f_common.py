"""Shared calibration and dual-axis (broken right axis) boxplot reader for Fig. S4F (PDF p. 18, 1200 dpi)."""

from __future__ import annotations

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (
    BrokenYCal,
    LinearCal,
    detect_ticks_along_spine,
    find_black_spines,
    find_horizontal_ticks_far,
    read_origin_boxplots,
    select_major_ticks,
)


def read_figs4f(image: np.ndarray, series: str) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    dark = gray < 105.0
    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]
    vr = v_spines[-1]

    # 3 categorical conditions: 1=VASP, 2=VASP noCP SLB, 3=VASP noCP ML
    x_cats = [1.0, 2.0, 3.0]
    r_plot_top = 40
    r_plot_bot = hr[0] - 2
    c_plot_left = vl[1] + 5
    c_plot_right = vr[0] - 5

    if series == "width":
        # Left y-axis: 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0 um (7 major ticks)
        yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=48, min_len=12, coord_min=40, coord_max=hr[0] - 5)
        yt = select_major_ticks(yt_raw, expected=7, frac=0.72)
        y_cal = LinearCal(yt, [3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0], max_res_px=2.5)
        return read_origin_boxplots(
            image, x_cats, "red", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=True, min_box_area=120,
        )
    elif series == "pore":
        # Right y-axis (broken): lower 4, 3, 2, 1, 0 um^2; upper 19, 18.5 (minor), 18 um^2
        all_maj = find_horizontal_ticks_far(dark, vr, "right", far_start=18, far_end=38, r_min=40, r_max=hr[0] - 5, min_run=8)
        lower_ticks = all_maj[-5:]
        lower_cal = LinearCal(lower_ticks, [4.0, 3.0, 2.0, 1.0, 0.0], max_res_px=2.5)
        # Upper segment: tick 19 (~153.0) and minor tick 18.5 (~351.0) cleanly above the break slash (~550..590)
        r19 = all_maj[0]
        inner = detect_ticks_along_spine(dark, vr, "y", "right", min_offset=2, max_offset=45, min_len=10, coord_min=40, coord_max=int(r19 + 260))
        minor_185 = [pos for pos, _ in inner if r19 + 140 < pos < r19 + 250]
        if not minor_185:
            raise RuntimeError(f"Missing 18.5 minor tick on Fig S4F right axis: inner={inner}")
        upper_cal = LinearCal([r19, minor_185[0]], [19.0, 18.5], max_res_px=1.5)
        break_row = 0.5 * (upper_cal.to_px(18.0) + lower_cal.to_px(4.0))
        y_cal = BrokenYCal(lower_cal, upper_cal, break_row)
        return read_origin_boxplots(
            image, x_cats, "blue", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=True, min_box_area=120, break_gap_max=220,
        )
    else:
        raise ValueError(series)
