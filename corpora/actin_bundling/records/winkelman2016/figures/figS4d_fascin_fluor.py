"""Calibrated reader for Fig. S4D (TMR-fascin fluorescence intensity bar chart)."""

from __future__ import annotations

import numpy as np

import sys
from pathlib import Path

_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, h_lines, v_lines = C.find_l_spines(
        image, thresh=85, row_frac=0.30, col_frac=0.30, row_range=(20, image.shape[0] - 40), col_range=(120, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0 (x 10^4 a.u.)
    yt_raw = C.find_ticks(dark, y_spine, (100, int(x_spine[1]) + 15), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 7, tol=8.0, min_step=280.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [3.0e4, 2.5e4, 2.0e4, 1.5e4, 1.0e4, 0.5e4, 0.0])), max_residual_px=2.0)

    H, W = gray.shape
    c0 = int(y_spine[2]) + 15
    r_above = int(x_spine[1]) - 10
    r_below = int(x_spine[2]) + 3
    blk = (image[:, :, 0] < 55) & (image[:, :, 1] < 55) & (image[:, :, 2] < 55)

    # Bar 3 (Parallel Bundle) is a tall gray bar rising above x-spine;
    # Bar 2 (Single Filament) is a thin negative blue bar just below x-spine;
    # Bar 1 (Background) sits at y ~ 0 (centered symmetrically to the left of Bar 2)
    band_up = gray[r_above - 80 : r_above, c0 : W - 20] < 235
    runs_up = [r for r in C.runs_1d(band_up.mean(axis=0) > 0.7, gap=5) if (r[1] - r[0] + 1) >= 400]
    R, B = image[:, :, 0].astype(np.float32), image[:, :, 2].astype(np.float32)
    blue_dn = (B[r_below : r_below + 45, c0 : W - 20] > 100) & (B[r_below : r_below + 45, c0 : W - 20] - R[r_below : r_below + 45, c0 : W - 20] > 30)
    runs_dn = [r for r in C.runs_1d(blue_dn.any(axis=0), gap=40) if (r[1] - r[0] + 1) >= 400]
    if len(runs_up) != 1 or len(runs_dn) != 1:
        raise RuntimeError(f"Fig. S4D: expected 1 positive + 1 negative bar, got up={runs_up}, dn={runs_dn}")

    ca3, cb3 = runs_up[0][0] + c0, runs_up[0][1] + c0
    ca2, cb2 = runs_dn[0][0] + c0, runs_dn[0][1] + c0
    c_mid3 = 0.5 * (ca3 + cb3) + 0.5
    c_mid2 = 0.5 * (ca2 + cb2) + 0.5
    # Bar 1 (Background) sits one bar-pitch to the left of Bar 2
    pitch = c_mid3 - c_mid2
    c_mid1 = c_mid2 - pitch
    # Refine c_mid1 to exact vertical error bar stem near x_spine
    stem_cols = np.where(blk[int(x_spine[1]) - 50 : int(x_spine[1]) - 8, int(c_mid1) - 80 : int(c_mid1) + 80].mean(axis=0) > 0.4)[0]
    if len(stem_cols) > 0:
        c_mid1 = float(stem_cols.mean()) + int(c_mid1) - 80 + 0.5

    out = []

    # Bar 1: Background (bar height ~ 0 right on x-spine, error bar straddling 0)
    eb_mask1 = blk.copy()
    eb_mask1[int(x_spine[1]) - 2 : int(x_spine[2]) + 3, :] = False
    u1, l1 = C.find_vertical_error_bar(
        eb_mask1, c_mid1, x_spine[0], row_min=int(x_spine[0]) - 200, row_max=int(x_spine[0]) + 200,
        stem_half_w=12, max_gap=20, cap_half_w=35, min_cap_width=14, start_up_offset=4, start_dn_offset=4
    )
    bar1_row = 0.5 * (u1 + l1) if (u1 is not None and l1 is not None) else float(x_spine[0])
    out.append(
        C.make_point_reading(
            1.0, c_mid1, bar1_row, y_axis, upper_err_row=u1, lower_err_row=l1, mirror_lower=True,
            err_px=2.5 + y_axis.residual_px, note=f"Fig. S4D bar 1 (Background); cal_res_y={y_axis.residual_px:.2f}px"
        )
    )

    # Bar 2: Single Filament (slightly negative blue bar below x-spine)
    margin2 = int(0.12 * (cb2 - ca2))
    left2 = gray[int(x_spine[2]) + 1 : int(x_spine[2]) + 120, ca2 + margin2 : int(c_mid2) - 40]
    right2 = gray[int(x_spine[2]) + 1 : int(x_spine[2]) + 120, int(c_mid2) + 40 : cb2 - margin2]
    is_b2 = (left2.mean(axis=1) < 235) & (right2.mean(axis=1) < 235)
    bot_off = 0
    while bot_off + 1 < len(is_b2) and is_b2[bot_off + 1]:
        bot_off += 1
    row2 = float(int(x_spine[2]) + 1 + bot_off)
    u2, l2 = C.find_vertical_error_bar(
        eb_mask1, c_mid2, row2, row_min=int(x_spine[0]) - 150, row_max=int(x_spine[0]) + 150,
        stem_half_w=12, max_gap=25, cap_half_w=35, min_cap_width=14, start_up_offset=2, start_dn_offset=2
    )
    out.append(
        C.make_point_reading(
            2.0, c_mid2, row2, y_axis, upper_err_row=u2, lower_err_row=l2, mirror_lower=True,
            err_px=2.5 + y_axis.residual_px, note=f"Fig. S4D bar 2 (Single Filament, slightly negative); cal_res_y={y_axis.residual_px:.2f}px"
        )
    )

    # Bar 3: Parallel Bundle (tall gray bar)
    margin3 = int(0.12 * (cb3 - ca3))
    left3 = gray[:r_above, ca3 + margin3 : int(c_mid3) - 40]
    right3 = gray[:r_above, int(c_mid3) + 40 : cb3 - margin3]
    is_b3 = (left3.mean(axis=1) < 235) & (right3.mean(axis=1) < 235)
    top3 = r_above - 1
    while top3 > 0 and is_b3[top3 - 1]:
        top3 -= 1
    row3 = float(top3)
    err_up3 = gray < 100
    err_up3[int(row3) + 1 :, :] = False
    u3, _ = C.find_vertical_error_bar(
        err_up3, c_mid3, row3, row_min=20, row_max=int(row3), stem_half_w=12, max_gap=8, cap_half_w=35, min_cap_width=16, start_up_offset=2
    )
    err_dn3 = blk.copy()
    err_dn3[: int(row3), :] = False
    _, l3 = C.find_vertical_error_bar(
        err_dn3, c_mid3, row3, row_min=int(row3), row_max=r_above, stem_half_w=12, max_gap=8, cap_half_w=35, min_cap_width=16, start_dn_offset=2
    )
    out.append(
        C.make_point_reading(
            3.0, c_mid3, row3, y_axis, upper_err_row=u3, lower_err_row=l3, mirror_lower=True,
            err_px=2.5 + y_axis.residual_px, note=f"Fig. S4D bar 3 (Parallel Bundle, positive control); cal_res_y={y_axis.residual_px:.2f}px"
        )
    )
    return out
