"""Shared calibrated reader for Fig. 6C (100% stacked bar chart of Rods vs Asters abundance, mean +/- SD)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read_stacked(image: np.ndarray, which: str) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.35, col_frac=0.35)
    vl0 = vl[0]
    y_end = np.where(dark[:, int(round(vl0[0]))])[0].max()
    hb = min(hl, key=lambda h: abs(h[0] - y_end))

    yt_c = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 20), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_c, 6, tol=12.0)
    ya = C.verify_axis("linear", list(zip(yt, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0])), max_residual_px=5.0)

    r_top_lim = int(yt[0]) - 20
    r_bot = int(hb[1]) - 5
    c_left = int(vl0[2]) + 20

    rgb = image.astype(float)
    # Grey Rods fill (low chroma, gray ~190-220)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    rods_fill = (chroma < 20.0) & (gray > 175.0) & (gray < 232.0)
    err_dark = gray < 75.0

    # Detect the 5 bars near baseline (y in [2%, 15%])
    rb0, rb1 = int(ya.to_pixel(15.0)), int(ya.to_pixel(2.0))
    bar_cols = C.clusters(np.where(rods_fill[rb0:rb1, c_left:].mean(axis=0) > 0.6)[0] + c_left, gap=30)
    bar_cols = [b for b in bar_cols if (b[2] - b[1] + 1) >= 180]
    if len(bar_cols) != 5:
        raise RuntimeError(f"Expected 5 bars in Fig. 6C, found {bar_cols}")

    x_vals = [0.0, 50.0, 100.0, 200.0, 500.0]
    out = []
    for k, ((cx, bx0, bx1), xv) in enumerate(zip(bar_cols, x_vals)):
        bw = bx1 - bx0 + 1
        wl = slice(bx0 + 20, bx0 + int(0.32 * bw))
        wr = slice(bx1 - int(0.32 * bw), bx1 - 20)
        rows_rods = np.where(
            (rods_fill[r_top_lim:r_bot, wl].mean(axis=1) > 0.6)
            & (rods_fill[r_top_lim:r_bot, wr].mean(axis=1) > 0.6)
        )[0] + r_top_lim
        fill_top = int(rows_rods[0])
        border_rows = np.where(
            (err_dark[max(r_top_lim, fill_top - 22) : fill_top + 2, wl].mean(axis=1) > 0.6)
            & (err_dark[max(r_top_lim, fill_top - 22) : fill_top + 2, wr].mean(axis=1) > 0.6)
        )[0] + max(r_top_lim, fill_top - 22)
        top_row = float(border_rows.mean()) + 0.5 if len(border_rows) else float(fill_top)

        if k == 0:
            u_row, l_row = top_row, top_row
        else:
            # Find upper and lower horizontal caps around cx (width ~30..90 px, clear of bar outline)
            ec0, ec1 = int(round(cx - 0.25 * bw)), int(round(cx + 0.25 * bw))
            cap_rows = []
            for r in range(r_top_lim + 25, r_bot):
                if abs(r - top_row) <= 12:
                    continue
                rs = C.runs_1d(err_dark[r, ec0 : ec1 + 1], gap=3)
                if rs:
                    best = max(rs, key=lambda t: t[1] - t[0])
                    w = best[1] - best[0] + 1
                    if 28 <= w <= 95:
                        cap_rows.append(r)
            cap_cs = C.clusters(np.array(cap_rows, dtype=int), gap=4)
            u_cands = [c[0] for c in cap_cs if c[0] < top_row]
            l_cands = [c[0] for c in cap_cs if c[0] > top_row]
            if not u_cands or not l_cands:
                raise RuntimeError(f"Fig. 6C bar {k+1}: missing error cap (u={u_cands}, l={l_cands})")
            u_row, l_row = min(u_cands), max(l_cands)

        h_px = 2.5 + ya.residual_px
        rods_y = float(ya.to_value(top_row))
        rods_hi = float(ya.to_value(u_row))
        rods_lo = float(ya.to_value(l_row))
        rods_hi_out = float(ya.to_value(u_row - h_px))
        rods_lo_out = float(ya.to_value(l_row + h_px))
        y_a, y_b = float(ya.to_value(top_row - h_px)), float(ya.to_value(top_row + h_px))
        lo_m, hi_m = min(y_a, y_b), max(y_a, y_b)

        if which == "rods":
            out.append({
                "x": xv,
                "y": rods_y,
                "lo": lo_m,
                "hi": hi_m,
                "bar_lo": min(rods_lo, rods_hi),
                "bar_hi": max(rods_lo, rods_hi),
                "bar_lo_outer": min(rods_lo_out, rods_hi_out),
                "bar_hi_outer": max(rods_lo_out, rods_hi_out),
                "px": [float(cx), float(top_row)],
                "px_lo": float(min(u_row, l_row)),
                "px_hi": float(max(u_row, l_row)),
                "note": f"Fig6C Rods [Arp2/3]={xv:g}nM, ya_res={ya.residual_px:.2f}px",
            })
        elif which == "asters":
            # Asters (%) = 100 - Rods (%); upper green segment spans from top_row up to 100%
            ast_y = 100.0 - rods_y
            ast_lo = 100.0 - max(rods_lo, rods_hi)
            ast_hi = 100.0 - min(rods_lo, rods_hi)
            ast_lo_out = 100.0 - max(rods_lo_out, rods_hi_out)
            ast_hi_out = 100.0 - min(rods_lo_out, rods_hi_out)
            out.append({
                "x": xv,
                "y": ast_y,
                "lo": 100.0 - hi_m,
                "hi": 100.0 - lo_m,
                "bar_lo": min(ast_lo, ast_hi),
                "bar_hi": max(ast_lo, ast_hi),
                "bar_lo_outer": min(ast_lo_out, ast_hi_out),
                "bar_hi_outer": max(ast_lo_out, ast_hi_out),
                "px": [float(cx), float(top_row)],
                "px_lo": float(min(u_row, l_row)),
                "px_hi": float(max(u_row, l_row)),
                "note": f"Fig6C Asters (=100-Rods) [Arp2/3]={xv:g}nM; cross marks Rods/Asters interface, ya_res={ya.residual_px:.2f}px",
            })
        else:
            raise ValueError(f"Unknown series {which!r}")
    return out
