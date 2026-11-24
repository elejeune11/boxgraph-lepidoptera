"""Calibrated reader for Extended Data Fig. 2d (VASP construct partitioning bar chart, mean +/- s.e.m., n=3)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.35)
    hb, vl0 = hl[-1], vl[0]
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 15), "y", side=-1, band=(6, 30), frac=0.4, max_width=26)
    yt = C.match_ticks_regular(yt, 4, tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt, [12.0, 8.0, 4.0, 0.0])), max_residual_px=5.0)

    r_top_lim = int(yt[0]) - 20
    r_bot = int(hb[1]) - 3
    c_left = int(vl0[2]) + 25
    c_right = image.shape[1] - 20

    # Grey filled bars (RGB ~128,128,128 -> 95 < gray < 170) sitting right above the x-axis
    bar_fill = (gray > 95.0) & (gray < 175.0)
    # Detect the 4 bar horizontal spans near bottom of plot (y in [0.2, 1.2])
    r_base0, r_base1 = int(ya.to_pixel(1.2)), int(ya.to_pixel(0.2))
    bar_cols = C.clusters(np.where(bar_fill[r_base0:r_base1, c_left:c_right].mean(axis=0) > 0.7)[0] + c_left, gap=15)
    bar_cols = [b for b in bar_cols if (b[2] - b[1] + 1) >= 120]
    if len(bar_cols) != 4:
        raise RuntimeError(f"Expected 4 bars in ED Fig. 2d, found {bar_cols}")

    # Dark grey error-bar strokes (gray < 90, distinct from ~130 grey bar fill!)
    err_dark = gray < 90.0
    out = []
    for k, (cx, bx0, bx1) in enumerate(bar_cols):
        bw = bx1 - bx0 + 1
        # Outer wings of bar (away from centered error bar) give exact bar top row
        w_l = slice(bx0 + 15, bx0 + int(0.25 * bw))
        w_r = slice(bx1 - int(0.25 * bw), bx1 - 15)
        rows_bar = np.where((bar_fill[r_top_lim:r_bot, w_l].mean(axis=1) > 0.7) & (bar_fill[r_top_lim:r_bot, w_r].mean(axis=1) > 0.7))[0] + r_top_lim
        top_row = float(rows_bar.min())

        # Centered error bar cap above and below top_row (horizontal dark stroke width 35..120 px)
        ec0, ec1 = int(round(cx - 0.25 * bw)), int(round(cx + 0.25 * bw))
        cap_rows = []
        for r in range(r_top_lim, r_bot):
            rs_h = C.runs_1d(err_dark[r, ec0 : ec1 + 1], gap=3)
            if rs_h and max(b - a + 1 for a, b in rs_h) >= 35:
                cap_rows.append(r)
        cap_cs = C.clusters(np.array(cap_rows, dtype=int), gap=4)
        u_row = cap_cs[0][0] if len(cap_cs) >= 1 else None
        l_row = cap_cs[-1][0] if len(cap_cs) >= 2 else None
        out.append(
            C.make_bar_reading(
                x_val=float(k + 1),
                col_center=cx,
                top_row=top_row,
                y_axis=ya,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=False,
                err_px=2.5 + ya.residual_px,
                note=f"ED2d construct={k+1}, ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
