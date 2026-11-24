"""Shared calibration and boxplot reader for Fig. 3F (PDF p. 5, 1200 dpi)."""

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


def read_fig3f(image: np.ndarray, series: str) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    dark = gray < 105.0
    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.40, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]
    vr = v_spines[-1]

    # 5 categorical conditions: 1=VASP, 2=VASP_DeltaTetra, 3=alpha-actinin, 4=mDia1, 5=mDia1+alpha-actinin
    x_cats = [1.0, 2.0, 3.0, 4.0, 5.0]
    r_plot_top = 20
    r_plot_bot = hr[0] - 2
    c_plot_left = vl[1] + 5
    c_plot_right = vr[0] - 5

    if series == "width":
        # Left y-axis: 6, 5, 4, 3, 2, 1, 0 um
        yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=40, min_len=8, coord_min=20, coord_max=hr[0] - 10)
        yt = select_major_ticks(yt_raw, expected=7, frac=0.70)
        y_cal = LinearCal(yt, [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0], max_res_px=2.5)
        return read_origin_boxplots(
            image, x_cats, "red", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=True, min_box_area=120,
        )
    elif series == "pore":
        # Right y-axis (broken): lower 20,18,16,14,12,10,8,6,4,2,0 um^2; upper 46, 44 um^2
        all_maj = find_horizontal_ticks_far(dark, vr, "right", far_start=15, far_end=30, r_min=20, r_max=hr[0] - 5, min_run=6)
        lower_ticks = all_maj[-11:]
        lower_cal = LinearCal(lower_ticks, [20.0, 18.0, 16.0, 14.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0], max_res_px=2.5)
        # Upper segment: tick 46 (~97.5) and minor tick 45 (~160.5) / tick 44 (~223)
        # Also detect minor tick 45 from inner strip if needed, or far strip if tick 44 is clean
        upper_cands = [r for r in all_maj[:-11] if r < lower_ticks[0] - 90]
        if len(upper_cands) == 2:
            # Verify 44 wasn't pulled by break slash: check spacing vs lower scale (127.1 px per 2 units)
            r46, r44 = upper_cands[0], upper_cands[1]
            # Also check minor tick 45 at inner strip (5..13 px right of spine)
            inner_ticks = detect_ticks_along_spine(dark, vr, "y", "right", min_offset=2, max_offset=35, min_len=7, coord_min=20, coord_max=int(r46 + 100))
            minor_45 = [pos for pos, l in inner_ticks if r46 + 40 < pos < r46 + 90]
            if minor_45:
                upper_cal = LinearCal([r46, minor_45[0]], [46.0, 45.0], max_res_px=1.5)
            else:
                upper_cal = LinearCal([r46, r44], [46.0, 44.0], max_res_px=2.5)
        else:
            raise RuntimeError(f"Unexpected upper ticks on Fig 3F right axis: all_maj={all_maj}")
        break_row = 0.5 * (upper_cal.to_px(44.0) + lower_cal.to_px(20.0))
        y_cal = BrokenYCal(lower_cal, upper_cal, break_row)
        return read_origin_boxplots(
            image, x_cats, "blue", y_cal, r_plot_top, r_plot_bot, c_plot_left, c_plot_right,
            include_whiskers=True, min_box_area=120, break_gap_max=160,
        )
    else:
        raise ValueError(series)
