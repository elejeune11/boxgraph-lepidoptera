"""Calibrated reader for Fig. 6F (actin-normalized relative enrichment of VASP vs Arp2/3 in aster spokes)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.35, col_frac=0.35)
    vl0 = vl[0]
    y_end = np.where(dark[:, int(round(vl0[0]))])[0].max()
    hb = min(hl, key=lambda h: abs(h[0] - y_end))

    yt_c = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 20), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_c, 6, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=2.0)

    r_top_lim = int(yt[0]) - 10
    r_bot = int(hb[1]) - 5
    c_left = int(vl0[2]) + 20
    bar_fill = (gray > 175.0) & (gray < 228.0)
    err_dark = gray < 80.0

    # Detect the 2 grey bars near baseline (y in [0.01, 0.05])
    rb0, rb1 = int(ya.to_pixel(0.05)), int(ya.to_pixel(0.01))
    bar_cols = C.clusters(np.where(bar_fill[rb0:rb1, c_left:].mean(axis=0) > 0.6)[0] + c_left, gap=30)
    bar_cols = [b for b in bar_cols if (b[2] - b[1] + 1) >= 200]
    if len(bar_cols) != 2:
        raise RuntimeError(f"Expected 2 bars in Fig. 6F, found {bar_cols}")

    names = ["VASP", "Arp2/3"]
    out = []
    for k, ((cx, bx0, bx1), name) in enumerate(zip(bar_cols, names)):
        bw = bx1 - bx0 + 1
        wl = slice(bx0 + 20, bx0 + int(0.32 * bw))
        wr = slice(bx1 - int(0.32 * bw), bx1 - 20)
        rows_bar = np.where(
            (bar_fill[r_top_lim:r_bot, wl].mean(axis=1) > 0.6)
            & (bar_fill[r_top_lim:r_bot, wr].mean(axis=1) > 0.6)
        )[0] + r_top_lim
        fill_top = int(rows_bar[0])
        border_rows = np.where(
            (err_dark[max(r_top_lim, fill_top - 22) : fill_top + 2, wl].mean(axis=1) > 0.6)
            & (err_dark[max(r_top_lim, fill_top - 22) : fill_top + 2, wr].mean(axis=1) > 0.6)
        )[0] + max(r_top_lim, fill_top - 22)
        top_row = float(border_rows.mean()) + 0.5 if len(border_rows) else float(fill_top)

        # Trace vertical error-bar stem centered at cx (below the significance bracket at y ~ 0.82)
        r_search_top = int(ya.to_pixel(0.78))
        stem = err_dark[r_search_top:r_bot, int(round(cx)) - 12 : int(round(cx)) + 13].mean(axis=1) > 0.45
        s_idx = np.where(stem)[0] + r_search_top
        cs = C.clusters(s_idx, gap=6)
        seg = min(cs, key=lambda t: abs(0.5 * (t[1] + t[2]) - top_row))
        u_row = float(seg[1]) + 0.5
        l_row = float(seg[2]) + 0.5

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
                note=f"Fig6F {name} (x={k+1}), ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
