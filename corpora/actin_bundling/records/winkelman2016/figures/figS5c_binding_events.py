"""Calibrated reader for Fig. S5C (single-molecule Cy5-fascin binding event frequency bar chart)."""

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
        image, thresh=85, row_frac=0.30, col_frac=0.30, row_range=(20, image.shape[0] - 20), col_range=(120, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 4, 3, 2, 1, 0 (5 sits near row ~160 at top of y-spine)
    yt_raw = C.find_ticks(dark, y_spine, (400, int(x_spine[1]) + 15), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 5, tol=8.0, min_step=380.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [4.0, 3.0, 2.0, 1.0, 0.0])), max_residual_px=2.5)

    H, W = gray.shape
    R = image[:, :, 0].astype(np.float32)
    G = image[:, :, 1].astype(np.float32)
    B = image[:, :, 2].astype(np.float32)
    r_base = int(x_spine[1]) - 15
    c0 = int(y_spine[2]) + 15

    # Bar 2 (Fascin domain) is tall red bar
    red = (R > 140) & (R - G > 60) & (R - B > 60)
    red_cols = red[r_base - 200 : r_base, c0 : W - 20].mean(axis=0) > 0.5
    runs_r = [r for r in C.runs_1d(red_cols, gap=5) if (r[1] - r[0] + 1) >= 300]
    if len(runs_r) != 1:
        raise RuntimeError(f"Fig. S5C: expected 1 red Fascin domain bar, found {runs_r}")
    ca2, cb2 = runs_r[0][0] + c0, runs_r[0][1] + c0
    c_mid2 = 0.5 * (ca2 + cb2) + 0.5

    # Bar 1 (alpha-Actinin domain) sits to the left of Bar 2 with y ~ 0 (small dark tick/bar on x-spine)
    c_mid1 = 0.5 * (float(y_spine[2]) + float(ca2))
    # Refine c_mid1 onto the small vertical error bar / baseline bar at x_spine
    blk = (R < 55) & (G < 55) & (B < 55)
    stem1 = np.where(blk[int(x_spine[1]) - 50 : int(x_spine[1]) - 6, int(c_mid1) - 100 : int(c_mid1) + 100].mean(axis=0) > 0.3)[0]
    if len(stem1) > 0:
        c_mid1 = float(stem1.mean()) + int(c_mid1) - 100 + 0.5

    # Check if there is a tiny upper error bar above x_spine at c_mid1
    eb1_up = blk.copy()
    eb1_up[int(x_spine[1]) - 1 :, :] = False
    u1, _ = C.find_vertical_error_bar(
        eb1_up, c_mid1, float(x_spine[1]), row_min=int(x_spine[1]) - 120, row_max=int(x_spine[1]) - 2,
        stem_half_w=12, max_gap=6, cap_half_w=25, min_cap_width=10, start_up_offset=2
    )
    # Bar 1 height is visually ~0 (0.01-0.02, essentially zero events on alpha-actinin domains)
    bar1_row = float(x_spine[0])
    rd1 = C.make_point_reading(
        1.0, c_mid1, bar1_row, y_axis, upper_err_row=u1, lower_err_row=bar1_row, mirror_lower=False,
        err_px=2.5 + y_axis.residual_px, note=f"Fig. S5C bar 1 (alpha-Actinin domain, ~0 events); cal_res_y={y_axis.residual_px:.2f}px"
    )

    # Bar 2 top (scan upward inside red bar away from central black error bar)
    margin2 = int(0.12 * (cb2 - ca2))
    left_r = red[:r_base, ca2 + margin2 : int(c_mid2) - 40].mean(axis=1) > 0.6
    right_r = red[:r_base, int(c_mid2) + 40 : cb2 - margin2].mean(axis=1) > 0.6
    top2 = r_base - 1
    while top2 > 0 and (left_r[top2 - 1] and right_r[top2 - 1]):
        top2 -= 1
    row2 = float(top2)

    err_up2 = gray < 100
    err_up2[int(row2) + 1 :, :] = False
    u2, _ = C.find_vertical_error_bar(
        err_up2, c_mid2, row2, row_min=10, row_max=int(row2), stem_half_w=14, max_gap=8, cap_half_w=35, min_cap_width=14, start_up_offset=2
    )
    err_dn2 = blk.copy()
    err_dn2[: int(row2), :] = False
    _, l2 = C.find_vertical_error_bar(
        err_dn2, c_mid2, row2, row_min=int(row2), row_max=r_base, stem_half_w=14, max_gap=8, cap_half_w=35, min_cap_width=14, start_dn_offset=2
    )
    rd2 = C.make_point_reading(
        2.0, c_mid2, row2, y_axis, upper_err_row=u2, lower_err_row=l2, mirror_lower=True,
        err_px=2.5 + y_axis.residual_px, note=f"Fig. S5C bar 2 (Fascin domain); cal_res_y={y_axis.residual_px:.2f}px"
    )
    return [rd1, rd2]
