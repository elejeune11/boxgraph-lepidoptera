"""Shared calibrated reader for Fig. 4F (3 stacked bar panels of 3D morphology abundance: Shell, Ring, Disc)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read_panel(image: np.ndarray, row_idx: int, label: str) -> list[dict]:
    """Read one stacked subpanel (row_idx = 0 for 0 nM, 1 for 100 nM, 2 for 150 nM)."""
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.5, col_frac=0.2)
    if len(hl) != 6 or len(vl) < 2:
        raise RuntimeError(f"Expected 6 horizontal box borders and >=2 vertical spines in Fig. 4F, got h={hl}, v={vl}")
    r_top_box = hl[2 * row_idx]
    r_bot_box = hl[2 * row_idx + 1]
    vl_left = vl[0]
    vl_right = vl[-1]

    # Calibrate y-axis for this subpanel using its 5 ticks (100, 75, 50, 25, 0)
    yt_c = C.find_outward_ticks(
        dark, vl_left, (int(r_top_box[1]) - 15, int(r_bot_box[2]) + 15),
        "y", side=-1, band=(5, 28), frac=0.45, max_width=22,
    )
    yt = C.match_ticks_regular(yt_c, 5, tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt, [100.0, 75.0, 50.0, 25.0, 0.0])), max_residual_px=3.0)

    r_top_lim = int(r_top_box[2]) + 6
    r_bot = int(r_bot_box[1]) - 4
    c_left = int(vl_left[2]) + 30
    c_right = int(vl_right[1]) - 30

    # Grey bar fill (low chroma, gray 175..228); excludes green/pink replicate dots
    rgb = image.astype(float)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    bar_fill = (chroma < 20.0) & (gray > 175.0) & (gray < 228.0)
    err_dark = (chroma < 25.0) & (gray < 75.0)

    # Detect the 3 bars (Shell, Ring, Disc) near bottom of panel (y in [1%, 6%]); bridge colored replicate dots up to ~80 px wide
    rb0, rb1 = int(ya.to_pixel(6.0)), int(ya.to_pixel(1.0))
    bar_cols = C.clusters(np.where(bar_fill[rb0:rb1, c_left:c_right].mean(axis=0) > 0.35)[0] + c_left, gap=95)
    bar_cols = [b for b in bar_cols if (b[2] - b[1] + 1) >= 400]
    if len(bar_cols) != 3:
        raise RuntimeError(f"Fig. 4F ({label}): expected 3 bars, found {bar_cols}")

    cats = ["Shell", "Ring", "Disc"]
    out = []
    for k, ((cx, bx0, bx1), cat) in enumerate(zip(bar_cols, cats)):
        bw = bx1 - bx0 + 1
        wl = slice(bx0 + 25, bx0 + int(0.24 * bw))
        wr = slice(bx1 - int(0.24 * bw), bx1 - 25)
        rows_bar = np.where(
            (bar_fill[r_top_lim:r_bot, wl].mean(axis=1) > 0.55)
            & (bar_fill[r_top_lim:r_bot, wr].mean(axis=1) > 0.55)
        )[0] + r_top_lim
        fill_top = int(rows_bar[0])
        border_rows = np.where(
            (err_dark[max(r_top_lim, fill_top - 20) : fill_top + 2, wl].mean(axis=1) > 0.55)
            & (err_dark[max(r_top_lim, fill_top - 20) : fill_top + 2, wr].mean(axis=1) > 0.55)
        )[0] + max(r_top_lim, fill_top - 20)
        top_row = float(border_rows.mean()) + 0.5 if len(border_rows) else float(fill_top)

        # Upper black SD error-bar cap (thin horizontal stroke centered on cx; bridge colored dots up to ~90 px wide; ignore tall text labels)
        ec0, ec1 = int(round(cx - 0.34 * bw)), int(round(cx + 0.34 * bw))
        cap_rows = []
        for r in range(r_top_lim, int(top_row) - 12):
            rs = C.runs_1d(err_dark[r, ec0 : ec1 + 1], gap=95)
            if rs:
                best = max(rs, key=lambda t: t[1] - t[0])
                w = best[1] - best[0] + 1
                n_blk = int(err_dark[r, ec0 + best[0] : ec0 + best[1] + 1].sum())
                # Require solid cap stroke (high pixel density) so letter strokes in '100 nM Arp2/3' cannot pass
                if w >= 160 and n_blk >= 130 and (n_blk / w) >= 0.48:
                    cap_rows.append(r)
        cap_cs = [c for c in C.clusters(np.array(cap_rows, dtype=int), gap=4) if (c[2] - c[1] + 1) <= 22]
        if not cap_cs:
            raise RuntimeError(f"Fig. 4F ({label}) {cat}: upper SD cap not found")
        # Pick the cap cluster closest above top_row (below any panel header text)
        u_row = cap_cs[-1][0]
        out.append(
            C.make_bar_reading(
                x_val=float(k + 1),
                col_center=cx,
                top_row=top_row,
                y_axis=ya,
                upper_err_row=u_row,
                lower_err_row=None,
                mirror_lower=True,
                err_px=2.5 + ya.residual_px,
                note=f"Fig4F {label} {cat} (x={k+1}), ya_res={ya.residual_px:.2f}px",
            )
        )
    return out
