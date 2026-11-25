"""Calibrated reader for Fig. 5C (simulated CytoSim Relative actin enrichment vs [Arp2/3])."""

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
    ya = C.verify_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=2.5)

    # 8 evenly spaced x-ticks corresponding to 0, 10, 25, 50, 75, 100, 125, 150 nM
    xt_c = C.find_outward_ticks(dark, hb, (int(vl0[2]) + 20, image.shape[1] - 10), "x", side=+1, band=(5, 35), frac=0.45, max_width=24)
    xt = C.match_ticks_regular(xt_c, 8, tol=8.0)
    dx_bar = float(np.median(np.diff(xt)))

    r_top_lim = int(yt[0]) - 20
    r_bot = int(hb[1]) - 5
    bar_fill = (gray > 175.0) & (gray < 230.0)
    err_dark = gray < 80.0

    x_vals = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0]
    out = []
    for k, (cx, xv) in enumerate(zip(xt, x_vals)):
        if k == 0:
            # 0 nM Arp2/3 has zero RAE (no visible bar above baseline)
            out.append(
                C.make_bar_reading(
                    x_val=xv,
                    col_center=cx,
                    top_row=float(yt[-1]),
                    y_axis=ya,
                    upper_err_row=float(yt[-1]),
                    lower_err_row=float(yt[-1]),
                    mirror_lower=False,
                    err_px=3.0 + ya.residual_px,
                    note=f"Fig5C [Arp2/3]=0nM (zero height on baseline), ya_res={ya.residual_px:.2f}px",
                )
            )
            continue

        # Bar wings within +/- 0.38*dx_bar around x-tick center, clear of centered error bar
        hw = int(round(0.36 * dx_bar))
        bx0, bx1 = int(round(cx - hw)), int(round(cx + hw))
        wl = slice(bx0 + 14, int(round(cx - 0.16 * dx_bar)))
        wr = slice(int(round(cx + 0.16 * dx_bar)), bx1 - 14)
        rows_bar = np.where(
            (bar_fill[r_top_lim:r_bot, wl].mean(axis=1) > 0.5)
            & (bar_fill[r_top_lim:r_bot, wr].mean(axis=1) > 0.5)
        )[0] + r_top_lim
        if len(rows_bar) == 0:
            raise RuntimeError(f"Fig. 5C bar {k+1} ([Arp2/3]={xv}): bar fill not found")
        fill_top = int(rows_bar[0])
        border_rows = np.where(
            (err_dark[max(r_top_lim, fill_top - 22) : fill_top + 2, wl].mean(axis=1) > 0.5)
            & (err_dark[max(r_top_lim, fill_top - 22) : fill_top + 2, wr].mean(axis=1) > 0.5)
        )[0] + max(r_top_lim, fill_top - 22)
        top_row = float(border_rows.mean()) + 0.5 if len(border_rows) else float(fill_top)

        # Error bar caps (both upper and lower horizontal caps centered on cx, width ~30..100 px)
        ec0, ec1 = int(round(cx - 0.22 * dx_bar)), int(round(cx + 0.22 * dx_bar))
        cap_rows = []
        for r in range(r_top_lim, r_bot):
            if abs(r - top_row) <= 10:
                continue
            rs = C.runs_1d(err_dark[r, ec0 : ec1 + 1], gap=3)
            if rs:
                best = max(rs, key=lambda t: t[1] - t[0])
                w = best[1] - best[0] + 1
                if 32 <= w <= 115:
                    cap_rows.append(r)
        cap_cs = C.clusters(np.array(cap_rows, dtype=int), gap=4)
        u_cands = [c[0] for c in cap_cs if c[0] < top_row]
        l_cands = [c[0] for c in cap_cs if c[0] > top_row]
        if not u_cands or not l_cands:
            raise RuntimeError(f"Fig. 5C bar {k+1}: missing error bar cap (u={u_cands}, l={l_cands})")
        u_row = min(u_cands)
        l_row = max(l_cands)
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
                note=f"Fig5C [Arp2/3]={xv:g}nM, ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
