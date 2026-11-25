"""Calibrated reader for Fig. 4B (peripheral actin abundance % vs [Arp2/3], mean +/- SD, n=3)."""

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
    ya = C.verify_axis("linear", list(zip(yt, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0])), max_residual_px=3.0)

    r_top_lim = int(yt[0]) - 30
    r_bot = int(hb[1]) - 5
    c_left = int(vl0[2]) + 20
    bar_fill = (gray > 175.0) & (gray < 228.0)

    # Detect the 3 bars near the bottom (y in [2%, 15%])
    rb0, rb1 = int(ya.to_pixel(15.0)), int(ya.to_pixel(2.0))
    bar_cols = C.clusters(np.where(bar_fill[rb0:rb1, c_left:].mean(axis=0) > 0.6)[0] + c_left, gap=25)
    bar_cols = [b for b in bar_cols if (b[2] - b[1] + 1) >= 180]
    if len(bar_cols) != 3:
        raise RuntimeError(f"Expected 3 bars in Fig. 4B, found {bar_cols}")

    x_vals = [0.0, 100.0, 150.0]
    err_dark = gray < 85.0
    out = []
    for k, ((cx, bx0, bx1), xv) in enumerate(zip(bar_cols, x_vals)):
        bw = bx1 - bx0 + 1
        # Outer wings of bar (clear of centered dots and error bar stem)
        wl = slice(bx0 + 18, bx0 + int(0.24 * bw))
        wr = slice(bx1 - int(0.24 * bw), bx1 - 18)
        rows_bar = np.where(
            (bar_fill[r_top_lim:r_bot, wl].mean(axis=1) > 0.6)
            & (bar_fill[r_top_lim:r_bot, wr].mean(axis=1) > 0.6)
        )[0] + r_top_lim
        # Top of grey fill sits immediately below the black top outline (~7 px stroke -> center ~3.5 px above fill)
        fill_top = int(rows_bar[0])
        # Refine exact center of top border line across outer wings
        border_rows = np.where(
            (err_dark[max(r_top_lim, fill_top - 25) : fill_top + 2, wl].mean(axis=1) > 0.6)
            & (err_dark[max(r_top_lim, fill_top - 25) : fill_top + 2, wr].mean(axis=1) > 0.6)
        )[0] + max(r_top_lim, fill_top - 25)
        top_row = float(border_rows.mean()) + 0.5 if len(border_rows) else float(fill_top)

        # Upper horizontal SD error-bar cap (width >= 120 px, wider than replicate dots ~40 px)
        ec0, ec1 = int(round(cx - 0.35 * bw)), int(round(cx + 0.35 * bw))
        cap_rows = []
        for r in range(r_top_lim, int(top_row) - 10):
            rs = C.runs_1d(err_dark[r, ec0 : ec1 + 1], gap=4)
            if rs and max(b - a + 1 for a, b in rs) >= 120:
                cap_rows.append(r)
        cap_cs = C.clusters(np.array(cap_rows, dtype=int), gap=4)
        if not cap_cs:
            raise RuntimeError(f"Fig. 4B bar {k+1}: upper SD cap not found")
        u_row = cap_cs[0][0]
        out.append(
            C.make_bar_reading(
                x_val=xv,
                col_center=cx,
                top_row=top_row,
                y_axis=ya,
                upper_err_row=u_row,
                lower_err_row=None,
                mirror_lower=True,
                err_px=2.5 + ya.residual_px,
                note=f"Fig4B [Arp2/3]={xv:g}nM, ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
