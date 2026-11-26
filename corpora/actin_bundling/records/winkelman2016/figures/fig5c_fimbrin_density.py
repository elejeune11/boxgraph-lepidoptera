"""Calibrated reader for Fig. 5C (fimbrin molecules/um F-actin bar chart)."""

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
        image, thresh=85, row_frac=0.30, col_frac=0.30, row_range=(40, image.shape[0] - 80), col_range=(200, image.shape[1] - 20)
    )
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 4, 3, 2, 1 (0 sits near x-spine)
    yt_raw = C.find_ticks(dark, y_spine, (60, int(x_spine[1]) - 80), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 4, tol=8.0, min_step=300.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [4.0, 3.0, 2.0, 1.0])), max_residual_px=2.0)

    H, W = gray.shape
    c0 = int(y_spine[2]) + 15
    r_base = int(x_spine[1]) - 15
    band = gray[r_base - 80 : r_base, c0 : W - 20] < 235
    runs = [r for r in C.runs_1d(band.mean(axis=0) > 0.7, gap=3) if (r[1] - r[0] + 1) >= 250]
    if len(runs) != 3:
        raise RuntimeError(f"Fig. 5C: expected 3 bars, found {runs}")

    labels = ["Single filaments (n=6)", "Parallel alpha-Actinin (n=23)", "Parallel Fascin (n=16)"]
    out = []
    for idx, (a, b) in enumerate(runs, start=1):
        ca, cb = a + c0, b + c0
        c_mid = 0.5 * (ca + cb) + 0.5
        margin = int(0.12 * (cb - ca))
        left_strip = gray[:r_base, ca + margin : int(c_mid) - 35]
        right_strip = gray[:r_base, int(c_mid) + 35 : cb - margin]
        is_bar = (left_strip.mean(axis=1) < 235) & (right_strip.mean(axis=1) < 235)
        top_r = r_base - 1
        while top_r > 0 and is_bar[top_r - 1]:
            top_r -= 1
        top_row = float(top_r)

        err_up = gray < 100
        err_up[int(top_row) + 1 :, :] = False
        u_row, _ = C.find_vertical_error_bar(
            err_up, c_mid, top_row, row_min=10, row_max=int(top_row), stem_half_w=12, max_gap=6, cap_half_w=35, min_cap_width=18, start_up_offset=2
        )
        # Error bar inside colored bar is pure black (R,G,B < 50)
        err_dn = (image[:, :, 0] < 55) & (image[:, :, 1] < 55) & (image[:, :, 2] < 55)
        err_dn[: int(top_row), :] = False
        _, l_row = C.find_vertical_error_bar(
            err_dn, c_mid, top_row, row_min=int(top_row), row_max=r_base, stem_half_w=12, max_gap=6, cap_half_w=35, min_cap_width=18, start_dn_offset=2
        )
        out.append(
            C.make_point_reading(
                float(idx),
                c_mid,
                top_row,
                y_axis,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=True,
                err_px=2.0 + y_axis.residual_px,
                note=f"Fig. 5C bar {idx} ({labels[idx - 1]}); cal_res_y={y_axis.residual_px:.2f}px",
            )
        )
    return out
