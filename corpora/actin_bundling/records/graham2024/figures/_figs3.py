"""Shared calibrated reader for Supplemental Fig. S3 (grouped bar chart of Fig. 4F with full SD bars)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C

CAT_RGB = [
    ("shell", np.array([176.0, 216.0, 168.0])),
    ("ring", np.array([208.0, 160.0, 200.0])),
    ("disc", np.array([200.0, 112.0, 168.0])),
]


def read_group(image: np.ndarray, group_idx: int, label: str) -> list[dict]:
    """Read one [Arp2/3] group (0=0nM, 1=100nM, 2=150nM) -> Shell (1), Ring (2), Disc (3)."""
    dark, gray, hl, vl = C.find_l_spines(image, thresh=100.0, row_frac=0.35, col_frac=0.35)
    vl0 = vl[0]
    y_end = np.where(dark[:, int(round(vl0[0]))])[0].max()
    hb = min(hl, key=lambda h: abs(h[0] - y_end))

    yt_c = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 20), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_c, 6, tol=14.0)
    ya = C.verify_axis("linear", list(zip(yt, [100.0, 80.0, 60.0, 40.0, 20.0, 0.0])), max_residual_px=6.0)

    r_top_lim = int(yt[0]) - 15
    r_bot = int(hb[1]) - 4
    c_left = int(vl0[2]) + 15
    rgb = image.astype(float)
    err_dark = gray < 95.0  # dark grey error bar strokes (RGB ~ 60-80)

    # Locate the 3 bars of each color across the 3 concentration groups near baseline (y in [1%, 6%])
    rb0, rb1 = int(ya.to_pixel(6.0)), int(ya.to_pixel(1.0))
    out = []
    for cat_i, (cat_name, target_rgb) in enumerate(CAT_RGB):
        mask = np.linalg.norm(rgb - target_rgb[None, None, :], axis=2) < 30.0
        cols = C.clusters(np.where(mask[rb0:rb1, c_left:].mean(axis=0) > 0.35)[0] + c_left, gap=35)
        cols = [b for b in cols if (b[2] - b[1] + 1) >= 120]
        if len(cols) != 3:
            raise RuntimeError(f"Fig. S3 {cat_name}: expected 3 bars across groups, found {cols}")
        cx, bx0, bx1 = cols[group_idx]
        bw = bx1 - bx0 + 1
        wl = slice(bx0 + 14, bx0 + int(0.32 * bw))
        wr = slice(bx1 - int(0.32 * bw), bx1 - 14)
        rows_bar = np.where(
            (mask[r_top_lim:r_bot, wl].mean(axis=1) > 0.5)
            & (mask[r_top_lim:r_bot, wr].mean(axis=1) > 0.5)
        )[0] + r_top_lim
        fill_top = int(rows_bar[0])
        border_rows = np.where(
            (err_dark[max(r_top_lim, fill_top - 18) : fill_top + 2, wl].mean(axis=1) > 0.5)
            & (err_dark[max(r_top_lim, fill_top - 18) : fill_top + 2, wr].mean(axis=1) > 0.5)
        )[0] + max(r_top_lim, fill_top - 18)
        top_row = float(border_rows.mean()) + 0.5 if len(border_rows) else float(fill_top)

        # Trace vertical error-bar stem (and horizontal caps) centered on cx!
        # Exclude significance brackets above the stem by taking the vertical stem cluster spanning top_row
        c_stem = slice(int(round(cx)) - 10, int(round(cx)) + 11)
        stem_rows = np.where(err_dark[r_top_lim:r_bot, c_stem].mean(axis=1) > 0.35)[0] + r_top_lim
        cs = C.clusters(stem_rows, gap=5)
        seg = min(cs, key=lambda t: abs(0.5 * (t[1] + t[2]) - top_row))
        # Refine cap centers near top and bottom of seg (horizontal cap width 25..70 px)
        ec0, ec1 = int(round(cx - 0.25 * bw)), int(round(cx + 0.25 * bw))
        u_rows = [
            r for r in range(max(r_top_lim, seg[1] - 4), min(r_bot, seg[1] + 18))
            if (rs := C.runs_1d(err_dark[r, ec0 : ec1 + 1], gap=3)) and 22 <= max(b - a + 1 for a, b in rs) <= 85
        ]
        l_rows = [
            r for r in range(max(r_top_lim, seg[2] - 18), min(r_bot, seg[2] + 4))
            if (rs := C.runs_1d(err_dark[r, ec0 : ec1 + 1], gap=3)) and 22 <= max(b - a + 1 for a, b in rs) <= 85
        ]
        u_row = float(np.mean(u_rows)) + 0.5 if u_rows else float(seg[1]) + 0.5
        l_row = float(np.mean(l_rows)) + 0.5 if l_rows else float(seg[2]) + 0.5

        out.append(
            C.make_bar_reading(
                x_val=float(cat_i + 1),
                col_center=cx,
                top_row=top_row,
                y_axis=ya,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=False,
                err_px=2.5 + ya.residual_px,
                note=f"FigS3 {label} {cat_name} (x={cat_i+1}), ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
