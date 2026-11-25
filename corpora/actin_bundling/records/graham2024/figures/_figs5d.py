"""Shared calibrated reader for Supplemental Fig. S5d (CytoSim R.A.E. bar charts across [VCA]=0..10 uM)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read_panel(image: np.ndarray, vca_label: str) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.5, col_frac=0.45)
    if len(vl) < 2 or len(hl) < 2:
        raise RuntimeError(f"Fig. S5d ({vca_label}): expected full box spines, got hl={hl}, vl={vl}")
    vl_l, vl_r = vl[0], vl[-1]
    hb_t, hb_b = hl[0], hl[-1]

    # 6 outward y-ticks (1.0, 0.8, 0.6, 0.4, 0.2, 0.0) on left spine
    yt_c = C.find_outward_ticks(
        dark, vl_l, (int(hb_t[1]) - 15, int(hb_b[2]) + 15),
        "y", side=-1, band=(4, 28), frac=0.45, max_width=22,
    )
    yt = C.match_ticks_regular(yt_c, 6, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=4.0)

    # 8 evenly spaced x-ticks (0, 10, 25, 50, 75, 100, 125, 150 nM) on bottom spine
    xt_c = C.find_outward_ticks(
        dark, hb_b, (int(vl_l[2]) + 10, int(vl_r[1]) - 10),
        "x", side=+1, band=(4, 28), frac=0.45, max_width=22,
    )
    xt = C.match_ticks_regular(xt_c, 8, tol=8.0)
    dx_bar = float(np.median(np.diff(xt)))

    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    # Olive-green bars (G > R + 15 and G > B + 35)
    bar_grn = (G > R + 15.0) & (G > B + 35.0)
    # Red vertical error bars (R > G + 45 and R > B + 45)
    err_red = (R > G + 45.0) & (R > B + 45.0)

    r_top_lim = int(hb_t[2]) + 5
    r_bot = int(hb_b[1]) - 3
    x_vals = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0]

    out = []
    for k, (cx, xv) in enumerate(zip(xt, x_vals)):
        if k == 0:
            # 0 nM Arp2/3 has zero RAE on the baseline
            out.append(
                C.make_bar_reading(
                    x_val=xv,
                    col_center=cx,
                    top_row=float(yt[-1]),
                    y_axis=ya,
                    upper_err_row=float(yt[-1]),
                    lower_err_row=float(yt[-1]),
                    mirror_lower=False,
                    err_px=2.5 + ya.residual_px,
                    note=f"FigS5d {vca_label} [Arp2/3]=0nM (zero on baseline), ya_res={ya.residual_px:.2f}px",
                )
            )
            continue

        # Sample green bar wings clear of centered red error bar
        hw = int(round(0.32 * dx_bar))
        wl = slice(int(round(cx - hw)), int(round(cx - 0.12 * dx_bar)))
        wr = slice(int(round(cx + 0.12 * dx_bar)), int(round(cx + hw)))
        rows_bar = np.where(
            (bar_grn[r_top_lim:r_bot, wl].mean(axis=1) > 0.5)
            & (bar_grn[r_top_lim:r_bot, wr].mean(axis=1) > 0.5)
        )[0] + r_top_lim
        if len(rows_bar) == 0:
            raise RuntimeError(f"Fig. S5d ({vca_label}) bar {k+1}: green bar not found")
        top_row = float(rows_bar[0])

        # Red vertical error bar centered on cx
        c_err = slice(int(round(cx)) - 12, int(round(cx)) + 13)
        red_rows = np.where(err_red[r_top_lim:r_bot, c_err].any(axis=1))[0] + r_top_lim
        if len(red_rows) >= 3:
            u_row = float(red_rows[0]) + 0.5
            l_row = float(red_rows[-1]) + 0.5
        else:
            u_row, l_row = top_row, top_row

        out.append(
            C.make_bar_reading(
                x_val=xv,
                col_center=cx,
                top_row=top_row,
                y_axis=ya,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=False,
                err_px=2.5 + ya.residual_px,
                note=f"FigS5d {vca_label} [Arp2/3]={xv:g}nM, ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
