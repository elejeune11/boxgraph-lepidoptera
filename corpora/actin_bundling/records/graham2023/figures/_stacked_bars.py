"""Shared calibrated reader for stacked bar charts with replicate dots and upper s.d. error bars:
- Main Fig. 2d (3 boxes: 1:20, 1:10, 1:7.5; 4 bars: Uniform, AR<1.1, AR>1.1, Rod)
- Main Fig. 3d (3 boxes: 1:20, 1:10, 1:7.5; 3 bars: Shell, Ring, Disc)
- Main Fig. 6d (2 boxes: 3% PEG, 10% PEG; 4 bars: Uniform, AR<1.1, AR>1.1, Rod)
- Extended Data Fig. 4c (3 boxes: 1:20, 1:10, 1:7.5; 3 bars: Shell, Ring, Disc)
"""

from __future__ import annotations

import numpy as np
import _calib as C


def _detect_boxes(image: np.ndarray, n_boxes: int):
    dark, gray, hl, _ = C.find_l_spines(image, row_frac=0.5, col_frac=0.15)
    if len(hl) < 2 * n_boxes:
        raise RuntimeError(f"Expected {2 * n_boxes} horizontal box borders, found {hl}")
    # Keep the 2*n_boxes longest/main frame borders (sorted top-to-bottom)
    hl_sorted = sorted(hl, key=lambda x: x[0])[: 2 * n_boxes]
    boxes = []
    for i in range(n_boxes):
        top_sp = hl_sorted[2 * i]
        bot_sp = hl_sorted[2 * i + 1]
        r0 = int(top_sp[2]) + 1
        r1 = int(bot_sp[1]) - 1
        sub = dark[r0 : r1 + 1, :]
        vc = C.clusters(np.where(sub.mean(axis=0) > 0.75)[0], gap=4)
        if len(vc) < 2:
            raise RuntimeError(f"Could not find left/right spines for box {i}: {vc}")
        c_left, c_right = vc[0], vc[-1]
        yt = C.find_outward_ticks(dark, c_left, (r0 - 15, r1 + 15), "y", side=-1, band=(6, 30), frac=0.4, max_width=24)
        yt = C.match_ticks_regular(yt, 5, tol=5.0)
        ya = C.verify_axis("linear", list(zip(yt, [100.0, 75.0, 50.0, 25.0, 0.0])), max_residual_px=3.0)
        boxes.append((r0, r1, int(c_left[2]) + 1, int(c_right[1]) - 1, ya))
    return dark, gray, boxes


def _find_bar_columns(image: np.ndarray, boxes: list, n_bars: int) -> list[tuple[int, int]]:
    """Pool vertical bar side-walls across all stacked boxes to find exact (bx0, bx1) of each bar."""
    _, _, c0, c1, _ = boxes[0]
    W_in = c1 - c0 + 1
    slot_w = W_in / float(n_bars)
    rgb = image.astype(float)
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]) - np.minimum(
        np.minimum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]
    )
    # Neutral dark grey/charcoal bar outline (excluding colorful replicate dots)
    border = (gray < 120.0) & (sat < 25.0)

    cols = []
    for k in range(n_bars):
        sc0 = int(round(c0 + k * slot_w)) + 10
        sc1 = int(round(c0 + (k + 1) * slot_w)) - 10
        # Accumulate vertical border pixel counts across all stacked panels
        v_score = np.zeros(sc1 - sc0 + 1, dtype=float)
        for r0, r1, _, _, ya in boxes:
            # Restrict to lower 95% of box below top frame
            sub = border[r0 + 20 : r1 - 2, sc0 : sc1 + 1]
            v_score += sub.sum(axis=0)
        # Find strong vertical edge clusters (bar left and right walls); allow short bars (~12 px tall in Fig. 6d)
        thr = max(12.0, 0.30 * float(v_score.max()))
        edges = C.clusters(np.where(v_score >= thr)[0] + sc0, gap=6)
        # Filter to edges separated by realistic bar width (25%-85% of slot width)
        best_pair = None
        best_sum = -1.0
        for i in range(len(edges)):
            for j in range(i + 1, len(edges)):
                w = edges[j][0] - edges[i][0]
                if 0.25 * slot_w <= w <= 0.85 * slot_w:
                    s = v_score[int(edges[i][0]) - sc0] + v_score[int(edges[j][0]) - sc0]
                    if s > best_sum:
                        best_sum = s
                        best_pair = (int(edges[i][1]), int(edges[j][2]))
        if best_pair is None:
            raise RuntimeError(f"Could not detect bar walls for bar {k} in slot [{sc0}, {sc1}]")
        cols.append(best_pair)
    return cols


def read_stacked_bars(
    image: np.ndarray,
    n_boxes: int,
    box_idx: int,
    n_bars: int,
    label: str,
) -> list[dict]:
    dark, gray, boxes = _detect_boxes(image, n_boxes)
    bar_cols = _find_bar_columns(image, boxes, n_bars)
    r0, r1, c0, c1, ya = boxes[box_idx]

    rgb = image.astype(float)
    sat = np.maximum(np.maximum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]) - np.minimum(
        np.minimum(rgb[:, :, 0], rgb[:, :, 1]), rgb[:, :, 2]
    )
    # Neutral dark grey/charcoal bar outline & error bars (never matches colored dots)
    border = (gray < 125.0) & (sat < 25.0)

    r_zero = float(ya.to_pixel(0.0))
    out = []
    for k, (bx0, bx1) in enumerate(bar_cols):
        bw = bx1 - bx0 + 1
        cx = 0.5 * (bx0 + bx1) + 0.5

        # Outer wings of the bar (immune to centered error bar and replicate dots)
        w_in = max(12, int(0.05 * bw))
        w_out = max(28, int(0.22 * bw))
        wing_l = slice(bx0 + w_in, bx0 + w_out)
        wing_r = slice(bx1 - w_out + 1, bx1 - w_in + 1)

        # Find the horizontal bar top border line across the bar wings (and full bar)
        top_row = None
        top_stroke_r0 = None
        for r in range(r0 + 8, int(round(r_zero)) + 2):
            wl_dark = border[r, wing_l].mean()
            wr_dark = border[r, wing_r].mean()
            full_dark = border[r, bx0 + 4 : bx1 - 3].mean()
            if (wl_dark > 0.65 and wr_dark > 0.65) or (full_dark > 0.72 and (wl_dark > 0.5 or wr_dark > 0.5)):
                r_end = r
                while r_end + 1 <= int(round(r_zero)) + 2 and (
                    (border[r_end + 1, wing_l].mean() > 0.5 and border[r_end + 1, wing_r].mean() > 0.5)
                    or border[r_end + 1, bx0 + 4 : bx1 - 3].mean() > 0.65
                ):
                    r_end += 1
                top_stroke_r0 = r
                top_row = 0.5 * (r + r_end) + 0.5
                break

        if top_row is None:
            top_row = r_zero
            top_stroke_r0 = int(round(r_zero))

        # Search strictly above top_stroke_r0 for the horizontal upper s.d. error cap.
        # A colored replicate dot may sit on the cap, splitting it into two dark segments;
        # checking outer tip span >= 0.34*bw and total dark pixels >= 0.16*bw catches both intact and dot-split caps
        # while rejecting single replicate dots (width ~40 px << 0.34*bw).
        cap_c0 = int(round(cx - 0.32 * bw))
        cap_c1 = int(round(cx + 0.32 * bw))
        err_top = None
        r_search_hi = max(r0 + 4, int(ya.to_pixel(102.0)))
        r_search_lo = top_stroke_r0 - 6
        if r_search_lo > r_search_hi:
            min_span = max(85, int(0.34 * bw))
            min_cnt = max(40, int(0.16 * bw))
            cap_rows = []
            for r in range(r_search_hi, r_search_lo + 1):
                dc = np.nonzero(border[r, cap_c0 : cap_c1 + 1])[0]
                if len(dc) >= min_cnt and (int(dc[-1] - dc[0] + 1) >= min_span):
                    cap_rows.append(r)
            if cap_rows:
                # Filter to thin horizontal cap lines (height <= 16 px, rejecting tall panel text)
                cap_cs = [c for c in C.clusters(np.array(cap_rows, dtype=int), gap=3) if (c[2] - c[1] + 1) <= 16]
                if cap_cs:
                    # Pick the cap closest above the bar top (avoids any stray upper text inside y<=102%)
                    err_top = cap_cs[-1][0]

        note = f"{label} cat={k+1}: bar=[{bx0},{bx1}], ya_res={ya.residual_px:.2f}px"
        if err_top is None:
            note += "; upper s.d. cap hidden inside/flush with bar top"
        rec = C.make_bar_reading(
            x_val=float(k + 1),
            col_center=cx,
            top_row=top_row,
            y_axis=ya,
            upper_err_row=err_top,
            lower_err_row=None,
            mirror_lower=True,
            err_px=2.0 + ya.residual_px,
            note=note,
        )
        # Clamp tiny negative zero-baseline jitter to 0.0
        if rec["y"] < 0.0 and rec["hi"] >= 0.0:
            rec["y"] = 0.0
            rec["lo"] = 0.0
        if "bar_lo" in rec and rec["bar_lo"] < 0.0:
            rec["bar_lo"] = max(0.0, rec["bar_lo"])
            rec["bar_lo_outer"] = max(0.0, rec["bar_lo_outer"])
        out.append(rec)
    return out
