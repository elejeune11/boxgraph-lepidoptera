"""Calibrated reader for Fig. S1A (leading and trailing barbed-end elongation rates, 5 bars)."""

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
        image, thresh=85, row_frac=0.30, col_frac=0.30, row_range=(450, 2150), col_range=(150, image.shape[1] - 20)
    )
    # x-spine sits near row ~2096 (avoiding p-value brackets above and underline below)
    x_spine = max(h_lines, key=lambda t: t[0])
    y_spine = min(v_lines, key=lambda t: t[0])

    # Inward y ticks: 10, 8, 6, 4, 2, 0 sub/s
    yt_raw = C.find_ticks(dark, y_spine, (450, int(x_spine[1]) + 15), "y", side=+1, band=(4, 24), frac=0.35, max_width=45)
    yt = C.match_ticks_regular(yt_raw, 6, tol=8.0, min_step=250.0)
    y_axis = C.verify_axis("linear", list(zip(yt, [10.0, 8.0, 6.0, 4.0, 2.0, 0.0])), max_residual_px=2.5)

    # Locate the 5 bars using their crisp black rectangular side borders right above x-spine!
    H, W = gray.shape
    c0 = int(y_spine[2]) + 15
    r_base = int(x_spine[1]) - 15
    # Scan columns in a mid-height band (y ~ 3..6 sub/s) where all 5 bars are present
    r_lo_b = int(round(y_axis.to_pixel(6.0)))
    r_hi_b = int(round(y_axis.to_pixel(2.0)))
    # Non-white pixels (including red/cyan stripes and black borders) across height
    nonwhite = (image[r_lo_b:r_hi_b, c0 : W - 20].min(axis=2) < 220).mean(axis=0) > 0.25
    runs = [r for r in C.runs_1d(nonwhite, gap=5) if (r[1] - r[0] + 1) >= 240]
    if len(runs) != 5:
        raise RuntimeError(f"Fig. S1A: expected 5 bars, found {runs}")

    labels = [
        "Leading barbed end (no error bar plotted)",
        "Trailing barbed end: Fascin alone",
        "Trailing barbed end: alpha-Actinin alone",
        "Trailing barbed end: Fascin domains (mixed)",
        "Trailing barbed end: alpha-Actinin domains (mixed)",
    ]
    blk = (image[:, :, 0] < 60) & (image[:, :, 1] < 60) & (image[:, :, 2] < 60)
    r_top_limit = int(round(y_axis.to_pixel(10.42)))  # strictly below the 'ns' significance bracket!

    out = []
    for idx, (a, b) in enumerate(runs, start=1):
        ca, cb = a + c0, b + c0
        c_mid = 0.5 * (ca + cb) + 0.5
        margin = int(0.10 * (cb - ca))
        # Every bar (solid or striped) has a solid horizontal black top border line!
        left_blk = blk[r_top_limit:r_base, ca + margin : int(c_mid) - 25].mean(axis=1) > 0.7
        right_blk = blk[r_top_limit:r_base, int(c_mid) + 25 : cb - margin].mean(axis=1) > 0.7
        top_cands = np.where(left_blk & right_blk)[0] + r_top_limit
        if len(top_cands) == 0:
            raise RuntimeError(f"Fig. S1A bar {idx}: no top black border found")
        top_cs = C.clusters(top_cands, gap=3)
        # For solid-black Bar 1, cluster spans entire bar height -> use top edge top_cs[0][1];
        # for outlined colored/striped bars 2-5, use center of the thin horizontal top border!
        top_row = float(top_cs[0][1]) + 2.5 if (top_cs[0][2] - top_cs[0][1]) > 30 else float(top_cs[0][0])

        u_row, l_row = None, None
        if idx > 1:
            # Trace black vertical SEM error bar above and below top_row (stop before significance bracket!)
            err_up = blk.copy()
            err_up[:r_top_limit, :] = False
            err_up[int(top_row) - 4 :, :] = False
            u_row, _ = C.find_vertical_error_bar(
                err_up, c_mid, top_row, row_min=r_top_limit, row_max=int(top_row) - 5, stem_half_w=10, max_gap=10, cap_half_w=28, min_cap_width=14, start_up_offset=6
            )
            err_dn = blk.copy()
            err_dn[: int(top_row) + 5, :] = False
            _, l_row = C.find_vertical_error_bar(
                err_dn, c_mid, top_row, row_min=int(top_row) + 5, row_max=r_base, stem_half_w=10, max_gap=10, cap_half_w=28, min_cap_width=14, start_dn_offset=6
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
                note=f"Fig. S1A bar {idx} ({labels[idx - 1]}); cal_res_y={y_axis.residual_px:.2f}px",
            )
        )
    return out
