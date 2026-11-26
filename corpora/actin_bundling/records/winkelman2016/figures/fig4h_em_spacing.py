"""Calibrated reader for Fig. 4H (negative-stain EM interfilament distance bar chart, 4 bars)."""

from __future__ import annotations

import numpy as np

import sys
from pathlib import Path

_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, h_lines, v_lines = C.find_l_spines(image, thresh=65, row_frac=0.35, col_frac=0.35)
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 30, 25, 20, 15, 10, 5 nm (0 sits at x-spine center)
    yt_raw = C.find_ticks(dark, y_spine, (250, int(x_spine[1]) - 30), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 6, tol=8.0, min_step=200.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [30.0, 25.0, 20.0, 15.0, 10.0, 5.0])), max_residual_px=2.5)

    H, W = gray.shape
    c0 = int(y_spine[2]) + 15
    r_base = int(x_spine[1]) - 15
    band = gray[r_base - 80 : r_base, c0:W] < 230
    runs = [r for r in C.runs_1d(band.mean(axis=0) > 0.7, gap=3) if (r[1] - r[0] + 1) >= 200]
    if len(runs) != 4:
        raise RuntimeError(f"Fig. 4H: expected 4 bars, found {runs}")

    labels = [
        "alpha-Actinin only",
        "alpha-Actinin (mixed with fascin)",
        "Fascin only",
        "Fascin (mixed with alpha-Actinin)",
    ]
    out = []
    for idx, (a, b) in enumerate(runs, start=1):
        ca, cb = a + c0, b + c0
        c_mid = 0.5 * (ca + cb) + 0.5
        margin = int(0.12 * (cb - ca))
        left_strip = gray[:r_base, ca + margin : int(c_mid) - 35]
        right_strip = gray[:r_base, int(c_mid) + 35 : cb - margin]
        is_bar = (left_strip.mean(axis=1) < 230) & (right_strip.mean(axis=1) < 230)
        top_r = r_base - 1
        while top_r > 0 and is_bar[top_r - 1]:
            top_r -= 1
        top_row = float(top_r)

        err_up = gray < 120
        err_up[int(top_row) + 1 :, :] = False
        u_row, _ = C.find_vertical_error_bar(
            err_up, c_mid, top_row, row_min=100, row_max=int(top_row), stem_half_w=10, max_gap=6, cap_half_w=25, min_cap_width=12, start_up_offset=2
        )
        bar_gray = float(np.median(left_strip[min(r_base - 10, int(top_row) + 25) : min(r_base, int(top_row) + 120)]))
        err_dn = gray < (0.55 * bar_gray)
        err_dn[: int(top_row), :] = False
        _, l_row = C.find_vertical_error_bar(
            err_dn, c_mid, top_row, row_min=int(top_row), row_max=r_base, stem_half_w=10, max_gap=6, cap_half_w=25, min_cap_width=12, start_dn_offset=2
        )
        if u_row is not None and (top_row - u_row) < 4:
            u_row = None
        if l_row is not None and (l_row - top_row) < 4:
            l_row = None
        # Mirror if bar 2 (very dark charcoal) or tiny bar hid the lower half
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
                note=f"Fig. 4H bar {idx} ({labels[idx - 1]}); cal_res_y={y_axis.residual_px:.2f}px",
            )
        )
    return out
