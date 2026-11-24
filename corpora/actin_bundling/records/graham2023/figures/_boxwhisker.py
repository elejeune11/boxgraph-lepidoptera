"""Shared calibrated reader for box-and-whisker aspect-ratio plots (Main Figs. 3g and 5a).

Reads the median (bisecting black line) as y, and the upper and lower s.d. whisker caps as
bar_hi and bar_lo, while logging Q1 and Q3 (box edges) in note.
"""

from __future__ import annotations

import numpy as np
import _calib as C


def read_boxwhisker(
    image: np.ndarray,
    y_ticks_vals: list[float],
    x_vals: list[float],
    y_max_search: float,
    label: str,
) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.35)
    hb = hl[-1]
    vl0 = vl[0]
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 15), "y", side=-1, band=(6, 30), frac=0.4, max_width=26)
    yt = C.match_ticks_regular(yt, len(y_ticks_vals), tol=5.0)
    ya = C.verify_axis("linear", list(zip(yt, y_ticks_vals)), max_residual_px=3.0)

    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    # Pure black box/whisker strokes (excludes teal individual data lines and peach box fill)
    black = (gray < 55.0) & (sat < 30.0)

    r_top_lim = max(0, int(ya.to_pixel(y_max_search)))
    r_bot_lim = int(hb[1]) - 4
    c_left = int(vl0[2]) + 20
    c_right = image.shape[1] - 20

    # Locate the n_cats horizontal box clusters from black pixels inside the plot
    sub_blk = black[r_top_lim : r_bot_lim + 1, c_left : c_right + 1]
    col_has = sub_blk.sum(axis=0) > 8
    cat_clusters = C.clusters(np.where(col_has)[0] + c_left, gap=80)
    # Filter to real box clusters (width >= 150 px)
    cat_clusters = [cc for cc in cat_clusters if (cc[2] - cc[1] + 1) >= 150]
    if len(cat_clusters) != len(x_vals):
        raise RuntimeError(f"Expected {len(x_vals)} box clusters in {label}, found {cat_clusters}")

    out = []
    for xv, (cx, bx0, bx1) in zip(x_vals, cat_clusters):
        bw = bx1 - bx0 + 1
        patch = black[r_top_lim : r_bot_lim + 1, bx0 : bx1 + 1]
        # Measure longest horizontal black run on each row
        row_span = np.zeros(patch.shape[0], dtype=int)
        for ri in range(patch.shape[0]):
            rs_h = C.runs_1d(patch[ri], gap=4)
            if rs_h:
                row_span[ri] = max(b_ - a_ + 1 for a_, b_ in rs_h)

        # Whisker caps: widest/narrower horizontal lines (>= 0.28 * bw)
        cap_rows = np.where(row_span >= int(0.28 * bw))[0] + r_top_lim
        full_rows = np.where(row_span >= int(0.82 * bw))[0] + r_top_lim
        if len(cap_rows) == 0 or len(full_rows) == 0:
            raise RuntimeError(f"Could not find box/whiskers for x={xv} in {label}")

        cap_cs = C.clusters(cap_rows, gap=3)
        full_cs = C.clusters(full_rows, gap=3)

        # Upper and lower whisker caps are the topmost and bottommost horizontal strokes
        if len(cap_cs) >= 2:
            r_hi = cap_cs[0][0]
            r_lo = cap_cs[-1][0]
        else:
            r_hi = float(cap_cs[0][1])
            r_lo = float(cap_cs[0][2])

        # Full-width horizontal strokes are Q3 (top), median (middle), Q1 (bottom)
        if len(full_cs) >= 3:
            r_q3 = full_cs[0][0]
            r_q1 = full_cs[-1][0]
            # Bisecting median line is the interior full-width stroke
            r_med = full_cs[len(full_cs) // 2][0]
            q_note = f"Q1={ya.to_value(r_q1):.3f}, Q3={ya.to_value(r_q3):.3f}"
        elif len(full_cs) == 2:
            # Narrow box where median merges with Q1 or Q3
            r_q3, r_q1 = full_cs[0][0], full_cs[-1][0]
            r_med = 0.5 * (r_q3 + r_q1)
            q_note = f"narrow box [{ya.to_value(r_q1):.3f}, {ya.to_value(r_q3):.3f}]"
        else:
            # Ultra-compact box (e.g. 1:20 AR ~ 1.05 in Fig. 5a where Q1/median/Q3 merge into one band)
            r_med = full_cs[0][0]
            q_note = f"compact box band {full_cs[0][1]}..{full_cs[0][2]}"

        h_px = 2.0 + ya.residual_px
        yv = float(ya.to_value(r_med))
        ya_, yb_ = float(ya.to_value(r_med - h_px)), float(ya.to_value(r_med + h_px))
        bar_hi = float(ya.to_value(r_hi))
        bar_lo = float(ya.to_value(r_lo))
        bar_hi_out = float(ya.to_value(r_hi - h_px))
        bar_lo_out = float(ya.to_value(r_lo + h_px))
        out.append(
            {
                "x": float(xv),
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "bar_lo": min(bar_lo, bar_hi),
                "bar_hi": max(bar_lo, bar_hi),
                "bar_lo_outer": min(bar_lo_out, bar_hi_out),
                "bar_hi_outer": max(bar_lo_out, bar_hi_out),
                "px": [float(cx), float(r_med)],
                "px_lo": float(min(r_hi, r_lo)),
                "px_hi": float(max(r_hi, r_lo)),
                "note": f"{label} x={xv:g}: median (whiskers=s.d.), {q_note}, ya_res={ya.residual_px:.2f}px",
            }
        )
    return out
