"""Shared calibrated reader for condensate aspect-ratio box/scatter plots
(Figs. 3E, 4C, 5I, 6J, 7F).

Each condition column has:
- Light blue-grey 25%-75% IQR box with dark grey outline
- Horizontal black median line spanning the box width
- Solid black square marker = overall Mean
- 3 cyan circles/diamonds = replicate means
- Light grey scatter diamonds = individual condensates
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _common as C


def read_ar_boxplot(
    img: np.ndarray,
    y_tick_vals: list[float],
    x_vals: list[float],
    metric: str,
    n_cols: int | None = None,
    select_indices: list[int] | None = None,
    y_tick_max_len: int = 75,
    max_res_px: float = 3.5,
    label: str = "",
) -> list[dict]:
    """Read `metric` ('mean' or 'median') from an aspect-ratio box/scatter plot.
    Also records Q25, Q75, replicate means, and max visible individual point in the note."""
    n_total = n_cols if n_cols is not None else len(x_vals)
    ax = C.find_l_axes(
        img,
        dark_thresh=115.0,
        min_h_frac=0.45,
        min_v_frac=0.5,
        max_x_tick_len=70,
        max_y_tick_len=y_tick_max_len,
    )
    xt = C.match_regular_ticks(ax["x_ticks"], n_total, tol=5.0)
    yt = C.match_regular_ticks(ax["y_ticks"], len(y_tick_vals), tol=5.0)
    ya = C.make_axis("linear", list(zip(yt, y_tick_vals)), max_res_px=max_res_px)

    rgb = img.astype(np.float64)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    gray = rgb.mean(axis=2)
    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)

    # Pure black ink (Mean square and Median horizontal line)
    pure_black = (gray < 45.0) & (sat < 20.0)
    # Cyan replicate mean markers: high G and B, very low R
    cyan = (G - R > 60.0) & (B - R > 60.0) & (G > 140.0) & (B > 140.0)
    # Blue-grey 25-75% box fill: slightly bluish/greenish grey (B >= R + 4, G >= R + 4, gray ~ 155-215)
    box_fill = (B - R >= 4.0) & (G - R >= 4.0) & (gray > 145.0) & (gray < 215.0) & (sat < 45.0)
    # Grey individual scatter diamonds (light grey ~175-225, low sat)
    any_scatter = (gray < 230.0) & (sat < 25.0)

    r_top_plot = max(0, int(round(ya.to_pixel(max(y_tick_vals)))) - 15)
    r_bot_plot = ax["bot"]["lo"] - 4
    spacing = (xt[-1] - xt[0]) / max(n_total - 1, 1) if n_total > 1 else (ax["bot"]["end"] - ax["left"]["hi"]) / 2.0
    half_w = int(round(0.32 * spacing))

    # Mask legend region in upper-left of plot (if present: legend text sits in top 32% of plot height,
    # left 55% of plot width — well clear of every Mean square in Figs. 3E, 4C, 5I, 6J, 7F).
    c_leg_r = ax["left"]["hi"] + int(round(0.55 * (ax["bot"]["end"] - ax["left"]["hi"])))
    r_leg_b = int(round(ya.to_pixel(min(y_tick_vals) + 0.68 * (max(y_tick_vals) - min(y_tick_vals)))))

    indices = select_indices if select_indices is not None else list(range(n_total))
    out = []
    for xv, ci in zip(x_vals, indices):
        cx = xt[ci]
        c0 = max(ax["left"]["hi"] + 5, int(round(cx - half_w)))
        c1 = min(img.shape[1] - 1, int(round(cx + half_w)))

        # Exclude upper-left legend text from pure_black & scatter in any column overlapping [left, c_leg_r]
        blk_col = pure_black[r_top_plot : r_bot_plot + 1, c0 : c1 + 1].copy()
        sc_col_full = any_scatter[r_top_plot : r_bot_plot + 1, c0 : c1 + 1].copy()
        if c0 < c_leg_r:
            leg_c_end = min(c1, c_leg_r) - c0 + 1
            leg_r_end = max(0, r_leg_b - r_top_plot)
            blk_col[:leg_r_end, :leg_c_end] = False
            sc_col_full[:leg_r_end, :leg_c_end] = False

        # 1. Detect the solid black Mean square (approx 24-75 px square centered right on cx).
        # Fill any enclosed holes (e.g. a cyan replicate dot sitting inside the black square in Fig. 3E 0.5 uM),
        # then apply separable 8x8 opening (erases thin horizontal median lines <= 6 px while preserving small squares).
        sub_blk_filled = ndimage.binary_fill_holes(blk_col)
        er = ndimage.binary_erosion(sub_blk_filled, structure=np.ones((9, 1), dtype=bool))
        er = ndimage.binary_erosion(er, structure=np.ones((1, 9), dtype=bool))
        sq_mask = ndimage.binary_dilation(er, structure=np.ones((9, 1), dtype=bool))
        sq_mask = ndimage.binary_dilation(sq_mask, structure=np.ones((1, 9), dtype=bool))
        lab_sq, n_sq = ndimage.label(sq_mask)
        mean_row = None
        mean_c = cx
        mean_half = 2.0
        best_sq_score = None
        for li in range(1, n_sq + 1):
            rr_i, cc_i = np.nonzero(lab_sq == li)
            area = len(rr_i)
            h_sq = int(rr_i.max() - rr_i.min() + 1)
            w_sq = int(cc_i.max() - cc_i.min() + 1)
            c_cent = cc_i.mean() + c0 + 0.5
            if abs(c_cent - cx) <= 0.14 * spacing and area >= 80 and w_sq <= 95 and h_sq <= 95:
                # Re-measure top & bottom edges from unoccluded pure-black columns of the square
                sq_c0, sq_c1 = int(cc_i.min()), int(cc_i.max())
                sq_r0 = max(0, int(rr_i.min()) - 45)
                sq_r1 = min(sub_blk_filled.shape[0] - 1, int(rr_i.max()) + 45)
                col_spans = []
                for cc_k in range(sq_c0, sq_c1 + 1):
                    ys_k = np.nonzero(sub_blk_filled[sq_r0 : sq_r1 + 1, cc_k])[0] + sq_r0
                    if len(ys_k) == 0:
                        continue
                    for cen_k, s_k, e_k in C.clusters_1d(ys_k, gap=2):
                        if e_k >= int(rr_i.min()) and s_k <= int(rr_i.max()):
                            col_spans.append((e_k - s_k + 1, 0.5 * (s_k + e_k) + 0.5))
                if col_spans:
                    max_h = max(h for h, _ in col_spans)
                    tall_centers = [cen for h, cen in col_spans if h >= max_h - 4]
                    r_cent = float(np.median(tall_centers)) + r_top_plot
                else:
                    r_cent = rr_i.mean() + r_top_plot + 0.5
                score = abs(c_cent - cx) - 0.02 * area
                if best_sq_score is None or score < best_sq_score:
                    best_sq_score = score
                    mean_row = r_cent
                    mean_c = c_cent
                    mean_half = 1.5

        # Fallback if a cyan replicate dot covers the center of a smaller Mean square (leaving only its black perimeter):
        # combine blk_col with cyan strictly within |c - cx| <= 28, open with 9x9, and measure the pure-black vertical extent.
        if mean_row is None:
            cyn_col = cyan[r_top_plot : r_bot_plot + 1, c0 : c1 + 1].copy()
            if c0 < c_leg_r:
                cyn_col[: max(0, r_leg_b - r_top_plot), : min(c1, c_leg_r) - c0 + 1] = False
            bc = blk_col.copy()
            c_rel_l = max(0, int(round(cx - 28)) - c0)
            c_rel_r = min(bc.shape[1] - 1, int(round(cx + 28)) - c0)
            bc[:, c_rel_l : c_rel_r + 1] |= cyn_col[:, c_rel_l : c_rel_r + 1]
            er2 = ndimage.binary_erosion(bc, structure=np.ones((9, 1), dtype=bool))
            er2 = ndimage.binary_erosion(er2, structure=np.ones((1, 9), dtype=bool))
            sq2 = ndimage.binary_dilation(er2, structure=np.ones((9, 1), dtype=bool))
            sq2 = ndimage.binary_dilation(sq2, structure=np.ones((1, 9), dtype=bool))
            lab2, n2 = ndimage.label(sq2)
            best_b = None
            for li in range(1, n2 + 1):
                rr_i, cc_i = np.nonzero(lab2 == li)
                c_cent = cc_i.mean() + c0 + 0.5
                if abs(c_cent - cx) > 0.14 * spacing or len(rr_i) < 200:
                    continue
                # Measure pure-black pixels belonging to this square (within +/- 6 px vertically, |c - cx| <= 28)
                r0_b = max(0, int(rr_i.min()) - 6)
                r1_b = min(blk_col.shape[0] - 1, int(rr_i.max()) + 6)
                b_rr, b_cc = np.nonzero(blk_col[r0_b : r1_b + 1, c_rel_l : c_rel_r + 1])
                if len(b_rr) >= 120:
                    r_cent = 0.5 * (float(b_rr.min()) + float(b_rr.max())) + r0_b + r_top_plot + 0.5
                    sc = -len(b_rr)
                    if best_b is None or sc < best_b:
                        best_b = sc
                        mean_row = r_cent
                        mean_c = cx
                        mean_half = 2.0

        r_col_top = r_leg_b if (ci == 0 and c0 < c_leg_r) else r_top_plot

        # 2. Detect the horizontal black Median line spanning the outer wings of the box
        # (ignoring the central +/- 55 px around cx where the Mean square & cyan dots sit!)
        box_half = int(round(0.24 * spacing))
        wl0 = max(c0, int(round(cx - box_half)))
        wl1 = max(wl0 + 5, int(round(cx - 58)))
        wr0 = min(c1 - 5, int(round(cx + 58)))
        wr1 = min(c1, int(round(cx + box_half)))

        dark_line = (gray < 85.0) & (sat < 22.0)
        med_rows = []
        for r in range(r_col_top, r_bot_plot + 1):
            lf = pure_black[r, wl0 : wl1 + 1].mean()
            rf = pure_black[r, wr0 : wr1 + 1].mean()
            if lf > 0.55 and rf > 0.55:
                med_rows.append((r, lf + rf))
        if not med_rows:
            for r in range(r_col_top, r_bot_plot + 1):
                lf = dark_line[r, wl0 : wl1 + 1].mean()
                rf = dark_line[r, wr0 : wr1 + 1].mean()
                if lf > 0.45 and rf > 0.45:
                    med_rows.append((r, lf + rf))
        if not med_rows:
            raise RuntimeError(f"{label} col={ci}: could not detect horizontal median line")

        med_clusters = C.clusters_1d(np.array([r for r, _ in med_rows]), gap=2)
        best_med = min(
            med_clusters,
            key=lambda cl: float(
                0.5
                * (
                    gray[cl[1] : cl[2] + 1, wl0 : wl1 + 1].mean()
                    + gray[cl[1] : cl[2] + 1, wr0 : wr1 + 1].mean()
                )
            ),
        )
        med_row = best_med[0]
        med_half = max(1.0, 0.5 * (best_med[2] - best_med[1] + 1))

        # 3. Measure 25%-75% IQR box extent from box_fill in the wings
        sub_bf = box_fill[r_col_top : r_bot_plot + 1, wl0 : wl1 + 1] | box_fill[r_col_top : r_bot_plot + 1, wr0 : wr1 + 1]
        bf_rows = np.flatnonzero(sub_bf.mean(axis=1) > 0.25) + r_col_top
        q_note = ""
        if len(bf_rows) >= 3:
            q75 = ya.to_value(float(bf_rows.min()))
            q25 = ya.to_value(float(bf_rows.max()))
            q_note = f"; IQR=[{q25:.3g},{q75:.3g}]"

        # 4. Measure cyan replicate mean dots
        sub_cy = cyan[r_col_top : r_bot_plot + 1, c0 : c1 + 1]
        lab_cy, n_cy = ndimage.label(sub_cy)
        cy_vals = []
        for li in range(1, n_cy + 1):
            rr_i, cc_i = np.nonzero(lab_cy == li)
            if len(rr_i) >= 40:
                cy_vals.append(ya.to_value(rr_i.mean() + r_col_top + 0.5))
        cy_vals.sort()
        cy_note = f"; rep_means=[{','.join(f'{v:.3g}' for v in cy_vals)}]" if cy_vals else ""

        # 5. Measure highest visible scatter point in this column (useful for Fig. 5I cutoff check!)
        sc_l = max(0, int(round(cx - 0.18 * spacing)) - c0)
        sc_r = min(sc_col_full.shape[1] - 1, int(round(cx + 0.18 * spacing)) - c0)
        sc_rows = np.flatnonzero(sc_col_full[:, sc_l : sc_r + 1].any(axis=1)) + r_top_plot
        max_note = f"; max_visible={ya.to_value(float(sc_rows.min())):.3g}" if len(sc_rows) > 0 else ""

        if metric == "mean":
            if mean_row is None:
                raise RuntimeError(f"{label} col={ci}: could not detect black Mean square")
            use_row, use_col, use_half = mean_row, mean_c, mean_half
            extra = f"median={ya.to_value(med_row):.3g}{q_note}{cy_note}{max_note}"
        elif metric == "median":
            use_row, use_col, use_half = med_row, cx, med_half
            m_str = f"mean={ya.to_value(mean_row):.3g}" if mean_row is not None else "mean=n/a"
            extra = f"{m_str}{q_note}{cy_note}{max_note}"
        else:
            raise ValueError(f"Unknown metric: {metric}")

        note = f"{label} {metric} col={ci+1}, ya_res={ya.residual_px:.2f}px; {extra}"
        out.append(
            C.reading_point(
                x_val=float(xv),
                col=use_col,
                row=use_row,
                y_axis=ya,
                half_px=use_half,
                note=note,
                clamp_min=0.0,
            )
        )
    return out
