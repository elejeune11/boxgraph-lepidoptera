"""Shared calibrated reader for box-and-whisker plots (Figs. 2B, 2D, 2F, 3C, 3E, S2).

Each box depicts IQR (Q1-Q3), median (bisecting line), mean ('x' or filled square),
and SD whiskers (symmetric caps at mean +/- SD). We read mean +/- SD (error_bars: sd)
and report Q1, median, and Q3 in each point's note.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read_boxplot(
    image: np.ndarray,
    y_ticks_vals: list[float],
    x_vals: list[float],
    label: str,
    max_y_for_boxes: float | None = None,
) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=100.0, row_frac=0.35, col_frac=0.35)
    # Pick the vertical y-axis spine (leftmost long vertical spine with >= 30% image height contiguous run)
    vl0 = None
    y_spine_top, y_spine_bot = 0, image.shape[0] - 1
    for v in vl:
        runs = C.runs_1d(dark[:, int(round(v[0]))], gap=3)
        if not runs:
            continue
        longest = max(runs, key=lambda t: t[1] - t[0])
        if longest[1] - longest[0] >= 0.35 * image.shape[0]:
            vl0 = v
            y_spine_top, y_spine_bot = longest
            break
    if vl0 is None:
        raise RuntimeError(f"{label}: could not locate y-axis spine")

    # Horizontal x-axis spine meets the bottom of the y-axis spine
    hb = min(hl, key=lambda h: abs(h[0] - y_spine_bot))

    yt_cand = C.find_outward_ticks(
        dark, vl0, (max(0, y_spine_top - 25), min(image.shape[0], int(hb[2]) + 25)),
        "y", side=-1, band=(5, 28), frac=0.45, max_width=24,
    )
    yt = C.match_ticks_regular(yt_cand, len(y_ticks_vals), tol=6.0)
    ya = C.verify_axis("linear", list(zip(yt, y_ticks_vals)), max_residual_px=3.0)

    # Pure black pixels (box frames, median lines, mean markers, SD whisker caps)
    rgb = image.astype(float)
    pure_blk = (rgb[:, :, 0] < 45.0) & (rgb[:, :, 1] < 45.0) & (rgb[:, :, 2] < 45.0)

    c_start = int(vl0[2]) + 25
    c_end = image.shape[1] - 10
    r_top = int(ya.to_pixel(max_y_for_boxes)) if max_y_for_boxes is not None else int(yt[0]) + 10
    r_bot = int(hb[1]) - 8

    # Detect box horizontal spans via wide pure-black lines (160..660 px wide)
    sub_b = pure_blk[r_top:r_bot, c_start:c_end]
    col_hits = np.zeros(sub_b.shape[1], dtype=int)
    for r in range(sub_b.shape[0]):
        for a, b in C.runs_1d(sub_b[r], gap=3):
            w = b - a + 1
            if 160 <= w <= 660:
                col_hits[a : b + 1] += 1
    b_cols = C.clusters(np.where(col_hits >= 5)[0] + c_start, gap=30)
    b_cols = [bc for bc in b_cols if (bc[2] - bc[1] + 1) >= 160]
    if len(b_cols) != len(x_vals):
        raise RuntimeError(f"{label}: expected {len(x_vals)} boxes, found {len(b_cols)}: {b_cols}")

    out = []
    for k, ((cx, x0, x1), xv) in enumerate(zip(b_cols, x_vals)):
        bw = x1 - x0 + 1
        wide_rows = []
        cap_rows = []
        for r in range(r_top, r_bot):
            row_seg = pure_blk[r, x0 : x1 + 1]
            rs = C.runs_1d(row_seg, gap=3)
            if not rs:
                continue
            best = max(rs, key=lambda t: t[1] - t[0])
            w = best[1] - best[0] + 1
            mid = 0.5 * (best[0] + best[1]) + x0
            if w >= 0.78 * bw:
                wide_rows.append(r)
            elif 0.25 * bw <= w <= 0.68 * bw and abs(mid - cx) < 0.15 * bw:
                cap_rows.append(r)

        wc = C.clusters(np.array(wide_rows, dtype=int), gap=3)
        cc = C.clusters(np.array(cap_rows, dtype=int), gap=3)
        if not wc:
            raise RuntimeError(f"{label} box {k+1}: no box horizontal lines found")

        q3_row = wc[0][0]
        q1_row = wc[-1][0]
        # Median line sits between Q3 and Q1 if 3 wide lines resolved; otherwise midpoint of tight box
        med_row = wc[1][0] if len(wc) >= 3 else 0.5 * (q3_row + q1_row)

        # Whisker caps: uppermost cap is at/above Q3, lowermost cap is at/below Q1
        upper_cands = [c[0] for c in cc if c[0] <= q3_row + 4.0]
        lower_cands = [c[0] for c in cc if c[0] >= q1_row - 4.0]
        u_row = min(upper_cands) if upper_cands else q3_row
        l_row = max(lower_cands) if lower_cands else q1_row

        # Mean position: midpoint of symmetric mean +/- SD whiskers; cross-check against interior black marker
        mean_row = 0.5 * (u_row + l_row)
        # Look for interior pure-black mean marker ('x' or square) when box is tall enough (> 25 px)
        if (q1_row - q3_row) > 25.0:
            c_mid0 = int(round(cx - 0.14 * bw))
            c_mid1 = int(round(cx + 0.14 * bw))
            # Search inside box excluding the horizontal Q3, median, Q1 lines (+/- 10 px)
            marker_rows = []
            for r in range(int(q3_row) + 10, int(q1_row) - 10):
                if any(abs(r - wline[0]) <= 10 for wline in wc):
                    continue
                cnt = int(pure_blk[r, c_mid0 : c_mid1 + 1].sum())
                if cnt >= 14:
                    marker_rows.append(r)
            if len(marker_rows) >= 12:
                mc = C.clusters(np.array(marker_rows, dtype=int), gap=6)
                # Pick cluster closest to symmetric whisker midpoint
                best_m = min(mc, key=lambda m: abs(m[0] - mean_row))
                if abs(best_m[0] - mean_row) <= 12.0:
                    mean_row = best_m[0]

        q1_val = ya.to_value(q1_row)
        med_val = ya.to_value(med_row)
        q3_val = ya.to_value(q3_row)
        sd_val = 0.5 * abs(ya.to_value(u_row) - ya.to_value(l_row))
        err_px = 2.0 + ya.residual_px
        out.append(
            C.make_bar_reading(
                x_val=float(xv),
                col_center=cx,
                top_row=mean_row,
                y_axis=ya,
                upper_err_row=u_row,
                lower_err_row=l_row,
                mirror_lower=False,
                err_px=err_px,
                note=(
                    f"{label} box={k+1}: mean+/-SD (SD={sd_val:.3f}), "
                    f"Q1={q1_val:.3f}, median={med_val:.3f}, Q3={q3_val:.3f}, ya_res={ya.residual_px:.2f}px"
                ),
            )
        )
    return out
