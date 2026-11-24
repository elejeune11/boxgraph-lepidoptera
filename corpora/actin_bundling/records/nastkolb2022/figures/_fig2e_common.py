"""Shared calibration and boxplot reader for Fig. 2E (PDF p. 4, 1200 dpi)."""

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


def read_fig2e(image: np.ndarray, series: str, include_whiskers: bool = True) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    H, W = gray.shape
    dark = gray < 105.0
    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]
    vr = v_spines[-1]

    # X axis: 0.00, 0.25, 0.50, 0.75, 1.00 uM
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=55, min_len=15, coord_min=vl[1], coord_max=vr[0])
    xt = select_major_ticks(xt_raw, expected=5, frac=0.72)
    x_cal = LinearCal(xt, [0.00, 0.25, 0.50, 0.75, 1.00], max_res_px=2.5)

    x_concs = [0.05, 0.125, 0.25, 0.50, 1.00]
    r_plot_top = 20
    r_plot_bot = hr[0] - 2
    c_plot_left = vl[1] + 5
    c_plot_right = vr[0] - 5

    if series == "width":
        # Left y-axis: 8, 7, 6, 5, 4, 3, 2, 1, 0 um
        yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=55, min_len=12, coord_min=20, coord_max=hr[0] - 10)
        yt = select_major_ticks(yt_raw, expected=9, frac=0.72)
        y_cal = LinearCal(yt, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0], max_res_px=2.5)
        return read_origin_boxplots(
            image, x_concs, "red", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=include_whiskers, min_box_area=120, x_cal=x_cal,
        )
    elif series == "pore":
        # Right y-axis (broken): lower 7..0 um^2, upper 18..16 um^2
        # Use far strip (25..44 px right of spine) to read ONLY major ticks and bypass minor ticks and break slash!
        all_maj = find_horizontal_ticks_far(dark, vr, "right", far_start=26, far_end=45, r_min=20, r_max=hr[0] - 5, min_run=8)
        # Filter ticks: upper segment is 18, 17, 16 (first 3) or 18, 17 if 16 overlaps break slash;
        # Lower segment is 7, 6, 5, 4, 3, 2, 1, 0 (last 8, evenly spaced ~228 px).
        lower_ticks = all_maj[-8:]
        lower_cal = LinearCal(lower_ticks, [7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0], max_res_px=2.5)
        upper_cands = [r for r in all_maj[:-8] if r < lower_ticks[0] - 400]
        if len(upper_cands) != 2:
            raise RuntimeError(f"Expected 2 clean upper major ticks (18, 17) on Fig 2E right axis, got {upper_cands} from {all_maj}")
        # Include minor tick 17.5 between 18 and 17 to verify linear residual on upper segment
        inner = detect_ticks_along_spine(dark, vr, "y", "right", min_offset=2, max_offset=45, min_len=12, coord_min=20, coord_max=int(upper_cands[1] + 20))
        minor_175 = [pos for pos, _ in inner if upper_cands[0] + 60 < pos < upper_cands[1] - 60]
        if minor_175:
            upper_cal = LinearCal([upper_cands[0], minor_175[0], upper_cands[1]], [18.0, 17.5, 17.0], max_res_px=2.0)
        else:
            upper_cal = LinearCal(upper_cands, [18.0, 17.0], max_res_px=2.0)
        # Break sits below tick 16 (~545) and above tick 7 (~810)
        break_row = 0.5 * (upper_cal.to_px(16.0) + lower_cal.to_px(7.0))
        y_cal = BrokenYCal(lower_cal, upper_cal, break_row)
        return read_origin_boxplots(
            image, x_concs, "blue", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=include_whiskers, min_box_area=120, break_gap_max=160, x_cal=x_cal,
        )
    else:
        raise ValueError(series)
