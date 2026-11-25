"""Shared calibrated stacked-bar reader for CytoSim bound-state fractions with staggered black SD stems
(Figs. 4K, 5D, S5E, S6C, S6D, S6G, S7E).

Handles both 1200-dpi vector main figures and 900-dpi rasterized Supplemental figures cleanly:
- Locates all n_rows stacked rectangular plot boxes automatically from dark frame spines.
- Calibrates each row's y-axis from its outward left y-ticks (or box inner top/bottom spines if ticks absent,
  as in Fig. S6D where y-tick labels are omitted and ticks align with box top=1.0, bottom=0.0).
- Separates stacked color segments per bar using median-RGB strips away from vertical black SD stems,
  splitting thin dark horizontal boundaries cleanly.
- Detects the horizontally staggered black SD vertical stem for each stacked segment (when present) to
  report exact bar_lo / bar_hi SD intervals!
"""

from __future__ import annotations

import numpy as np
import _common as C


def find_stacked_boxes(
    img: np.ndarray,
    n_rows: int,
    dark_thresh: float = 115.0,
    min_h_frac: float = 0.40,
    min_v_frac: float = 0.10,
    max_tick_len: int = 45,
) -> list[dict]:
    """Locate n_rows vertically stacked rectangular plot boxes, sorted top-to-bottom."""
    dark = img.sum(axis=2, dtype=np.uint16) < int(round(3.0 * dark_thresh))
    H, W = dark.shape
    hs = C.group_spines(dark, min_len=min_h_frac * W, gap=5, merge_dist=5)
    # Filter out any thin legend/text lines by requiring matched pairs sharing nearly identical (start, end)
    # All main plot box horizontal spines span the same column range [s_box, e_box]
    if len(hs) < 2 * n_rows:
        # Try slightly looser threshold
        hs = C.group_spines(dark, min_len=0.32 * W, gap=6, merge_dist=6)
    if len(hs) < 2 * n_rows:
        raise ValueError(f"Expected >= {2 * n_rows} horizontal box spines, found {len(hs)}: {hs}")

    # Group spines with similar horizontal extent (within 40 px of the longest spine)
    max_len = max(h["length"] for h in hs)
    box_hs = [h for h in hs if h["length"] >= 0.88 * max_len]
    # If panels touch/share a spine (len == n_rows + 1) vs separated (len == 2 * n_rows)
    box_hs.sort(key=lambda h: h["lo"])

    pairs: list[tuple[dict, dict]] = []
    if len(box_hs) == 2 * n_rows:
        for i in range(n_rows):
            pairs.append((box_hs[2 * i], box_hs[2 * i + 1]))
    elif len(box_hs) == n_rows + 1:
        for i in range(n_rows):
            pairs.append((box_hs[i], box_hs[i + 1]))
    else:
        # Pair consecutive spines whose gap > 60 px (reject legend bars above top box)
        # Take the bottom 2 * n_rows spines
        cand = box_hs[-2 * n_rows :]
        for i in range(n_rows):
            pairs.append((cand[2 * i], cand[2 * i + 1]))

    boxes = []
    for top_sp, bot_sp in pairs:
        r0 = top_sp["hi"] + 1
        r1 = bot_sp["lo"]
        # Find left & right vertical spines strictly within [top_sp['lo'], bot_sp['hi']]
        sub = dark[top_sp["lo"] : bot_sp["hi"] + 1, :]
        col_frac = sub.mean(axis=0)
        vc = C.clusters_1d(np.flatnonzero(col_frac > 0.72), gap=3)
        if len(vc) < 2:
            c_left_lo, c_left_hi = top_sp["start"], top_sp["start"] + 6
            c_right_lo, c_right_hi = top_sp["end"] - 6, top_sp["end"]
        else:
            c_left_lo, c_left_hi = vc[0][1], vc[0][2]
            c_right_lo, c_right_hi = vc[-1][1], vc[-1][2]
        c0 = c_left_hi + 1
        c1 = c_right_lo

        # Outward ticks along left y-spine and bottom x-spine
        y_strip = dark[:, max(0, c_left_lo - max_tick_len - 5) : c_left_lo][:, ::-1].T
        y_ticks = C.find_ticks_on_strip(
            y_strip, top_sp["lo"] - 6, bot_sp["hi"] + 7, min_len=3, max_len=max_tick_len, max_width=24
        )
        x_strip = dark[bot_sp["hi"] + 1 : min(H, bot_sp["hi"] + max_tick_len + 6), :]
        x_ticks = C.find_ticks_on_strip(
            x_strip, c_left_lo - 6, c_right_hi + 7, min_len=3, max_len=max_tick_len, max_width=24
        )
        boxes.append({
            "box": (int(c0), int(r0), int(c1), int(r1)),
            "top": top_sp,
            "bot": bot_sp,
            "left": (c_left_lo, c_left_hi, 0.5 * (c_left_lo + c_left_hi) + 0.5),
            "right": (c_right_lo, c_right_hi, 0.5 * (c_right_lo + c_right_hi) + 0.5),
            "y_ticks": y_ticks,
            "x_ticks": x_ticks,
        })
    return boxes


def _measure_staggered_sd(
    img: np.ndarray,
    bx0: int,
    bx1: int,
    r_top_box: int,
    r_bot_box: int,
    y_axis,
    cum_top_vals: list[float],
    seg_fracs: list[float],
    cls_idx: int,
    n_classes: int,
    stem_fracs: list[float] | None = None,
    black_thresh: float = 65.0,
) -> tuple[float | None, float | None, float]:
    """Locate the staggered vertical black SD bar for stacked class cls_idx inside bar [bx0, bx1].
    Returns (sd_val, stem_col, top_boundary_row_or_mid)."""
    bw = bx1 - bx0 + 1
    y_b = cum_top_vals[cls_idx]
    if stem_fracs is not None:
        frac_pos = stem_fracs[cls_idx]
    elif n_classes == 3:
        # Measured directly on Fig. 4K / S5E / S6D / S6G: left-shifted 0.1375 step from 0.225
        frac_pos = [0.225, 0.362, 0.500][cls_idx]
    elif n_classes == 5:
        # Measured directly on Fig. 5D / S6C: symmetric 0.1375 step centered at 0.500
        frac_pos = [0.225, 0.362, 0.500, 0.638, 0.775][cls_idx]
    elif n_classes == 6:
        # Measured directly on Fig. S7E: [0.207, 0.353, 0.494, 0.626, 0.772, 0.898]
        frac_pos = [0.208, 0.347, 0.486, 0.625, 0.764, 0.902][cls_idx]
    else:
        frac_pos = (cls_idx + 0.5) / n_classes

    target_rel_c = int(round(frac_pos * bw)) - 5
    # Absent segment (fraction < 0.0015) has no SD stem; return 0 cleanly
    if seg_fracs[cls_idx] < 0.0015:
        return 0.0, bx0 + frac_pos * bw, float(y_axis.to_pixel(y_b))

    sub = img[r_top_box + 2 : r_bot_box - 1, bx0 + 5 : bx1 - 4].astype(np.float64)
    if sub.shape[0] < 10 or sub.shape[1] < 10:
        return None, 0.5 * (bx0 + bx1), float(y_axis.to_pixel(y_b))

    gray = sub.mean(axis=2)
    sat = sub.max(axis=2) - sub.min(axis=2)
    blk = (gray < black_thresh) & (sat < 35.0)

    r_b_rel = int(round(y_axis.to_pixel(y_b))) - (r_top_box + 2)
    is_top_cls = (cls_idx == n_classes - 1)

    # Strictly < half the spacing between adjacent stems so a missing stem never grabs a neighbor's stem
    search_half = max(5, int(round(0.052 * bw)))
    sc0 = max(0, target_rel_c - search_half)
    sc1 = min(sub.shape[1] - 1, target_rel_c + search_half)

    # Exclude full-bar-width dark horizontal segment border lines BEFORE probing for the vertical stem!
    row_full_blk = blk.mean(axis=1) > 0.60
    blk_stem = blk & ~row_full_blk[:, None]

    # Look for a vertical black stem passing through r_b_rel (within +/- 10 px of y_b, or top 16 px for top class)
    if is_top_cls:
        r_probe0, r_probe1 = 0, min(sub.shape[0] - 1, 16)
    else:
        r_probe0 = max(0, r_b_rel - 10)
        r_probe1 = min(sub.shape[0] - 1, r_b_rel + 10)

    probe = blk_stem[r_probe0 : r_probe1 + 1, sc0 : sc1 + 1]
    col_hits = probe.sum(axis=0)
    hit_cols = np.flatnonzero(col_hits >= 2)
    if hit_cols.size == 0:
        return 0.0, bx0 + frac_pos * bw, float(y_axis.to_pixel(y_b))

    c_clusters = [c for c in C.clusters_1d(hit_cols, gap=2) if (c[2] - c[1] + 1) <= 24]
    if not c_clusters:
        return 0.0, bx0 + frac_pos * bw, float(y_axis.to_pixel(y_b))

    best_cl = min(c_clusters, key=lambda c: abs((sc0 + c[0]) - target_rel_c))
    c_l = max(0, sc0 + best_cl[1] - 1)
    c_r = min(sub.shape[1] - 1, sc0 + best_cl[2] + 1)
    stem_col = bx0 + 5 + sc0 + best_cl[0]

    stem_on = blk_stem[:, c_l : c_r + 1].any(axis=1)

    px_per_unit = abs(y_axis.to_pixel(0.0) - y_axis.to_pixel(1.0))

    if is_top_cls:
        # Walk DOWN from row 0 (y = 1.0) while stem continues
        dn = 0
        miss = 0
        for rr in range(0, sub.shape[0]):
            if stem_on[rr]:
                dn = rr
                miss = 0
            else:
                miss += 1
                if miss > 4:
                    break
        if dn <= 2:
            return 0.0, stem_col, float(y_axis.to_pixel(y_b))
        # Downward distance from y = 1.0 gives SD directly!
        r_bot_abs = (r_top_box + 2) + dn + 0.5
        sd_val = max(0.0, 1.0 - y_axis.to_value(r_bot_abs))
        return sd_val, stem_col, r_bot_abs
    else:
        # Walk UP and DOWN from r_b_rel
        r_seed = max(0, min(sub.shape[0] - 1, r_b_rel))
        up = r_seed
        miss = 0
        for rr in range(r_seed, -1, -1):
            if stem_on[rr] or row_full_blk[rr]:
                if stem_on[rr]:
                    up = rr
                miss = 0
            else:
                miss += 1
                if miss > 4:
                    break
        dn = r_seed
        miss = 0
        for rr in range(r_seed, sub.shape[0]):
            if stem_on[rr] or row_full_blk[rr]:
                if stem_on[rr]:
                    dn = rr
                miss = 0
            else:
                miss += 1
                if miss > 4:
                    break
        span_px = max(0, dn - up)
        if span_px <= 3:
            return 0.0, stem_col, float(y_axis.to_pixel(y_b))
        # Check if stem clipped at top (y=1.0) or bottom (y=0.0)
        clipped_top = up <= 2
        clipped_bot = dn >= sub.shape[0] - 3
        if clipped_top and not clipped_bot:
            sd_px = max(0.0, dn - r_b_rel)
        elif clipped_bot and not clipped_top:
            sd_px = max(0.0, r_b_rel - up)
        else:
            sd_px = 0.5 * span_px
        sd_val = sd_px / px_per_unit
        return sd_val, stem_col, float(y_axis.to_pixel(y_b))


def read_stacked_bars_panel(
    img: np.ndarray,
    n_rows: int,
    row_idx: int,
    x_vals: list[float],
    palette_rgb: list[tuple[float, float, float]],
    cls_idx: int,
    y_tick_vals: list[float] = (1.0, 0.5, 0.0),
    active_bar_indices: list[int] | None = None,
    n_slots: int | None = None,
    color_tol: float = 48.0,
    max_res_px: float = 4.0,
    has_sd: bool = True,
    use_box_y_fallback: bool = False,
    stem_fracs: list[float] | None = None,
    label: str = "",
) -> list[dict]:
    """Read one stacked class `cls_idx` across all bars in row `row_idx` of an `n_rows` stacked bar grid."""
    boxes = find_stacked_boxes(img, n_rows)
    b = boxes[row_idx]
    c0, r0, c1, r1 = b["box"]

    # Calibrate y-axis on left outward ticks (with fallback to box top/bottom spines if labels/ticks omitted)
    if not use_box_y_fallback and len(b["y_ticks"]) >= len(y_tick_vals):
        yt = C.match_regular_ticks(b["y_ticks"], len(y_tick_vals), tol=6.0)
        ya = C.make_axis("linear", list(zip(yt, y_tick_vals)), max_res_px=max_res_px)
    else:
        ya = C.make_axis("linear", [(float(r0), 1.0), (float(r1), 0.0)], max_res_px=max_res_px)

    # Bar centers from bottom panel's x_ticks (shared vertically across aligned rows)
    n_total = n_slots if n_slots is not None else len(x_vals)
    xt_all = boxes[-1]["x_ticks"]
    if len(xt_all) >= n_total:
        x_centers = C.match_regular_ticks(xt_all, n_total, tol=6.0)
    else:
        # Evenly divide box width if x-ticks are absent
        step = (c1 - c0) / float(n_total)
        x_centers = [c0 + (k + 0.5) * step for k in range(n_total)]

    bar_indices = active_bar_indices if active_bar_indices is not None else list(range(n_total))
    pal = np.asarray(palette_rgb, dtype=np.float64)
    n_cls = len(palette_rgb)
    row_0 = float(ya.to_pixel(0.0))
    row_1 = float(ya.to_pixel(1.0))
    px_per_unit = abs(row_0 - row_1)
    # Bar width estimation: spacing between centers * 0.72
    if n_total >= 2:
        bar_sp = (x_centers[-1] - x_centers[0]) / float(n_total - 1)
    else:
        bar_sp = float(c1 - c0)
    bar_half_w = int(round(0.36 * bar_sp))

    out = []
    for xv, bi in zip(x_vals, bar_indices):
        cx = x_centers[bi]
        bx0 = max(c0 + 2, int(round(cx - bar_half_w)))
        bx1 = min(c1 - 2, int(round(cx + bar_half_w)))
        # Refine exact bar left/right fill edges inside [bx0, bx1]
        sub_bar = img[r0 + 2 : r1 - 1, bx0 : bx1 + 1].astype(np.float64)
        # Classify every pixel in bar strip to nearest palette color, rejecting white background and black stems/borders
        gray = sub_bar.mean(axis=2)
        sat = sub_bar.max(axis=2) - sub_bar.min(axis=2)
        is_black = (gray < 45.0) & (sat < 22.0)
        is_white = (gray > 238.0) & (sat < 15.0)
        valid_px = ~is_black & ~is_white

        # Distances to palette
        dists = np.sqrt(((sub_bar[:, :, None, :] - pal[None, None, :, :]) ** 2).sum(axis=3))
        nearest = np.argmin(dists, axis=2)
        min_d = np.min(dists, axis=2)
        matched = valid_px & (min_d <= color_tol)

        # Tighten [bx0, bx1] to actual colored bar columns (where >= 30% of rows match any palette color)
        col_has_bar = matched.mean(axis=0) > 0.30
        bar_cols_rel = np.flatnonzero(col_has_bar)
        if bar_cols_rel.size >= 10:
            bx0_ref = bx0 + int(bar_cols_rel[0])
            bx1_ref = bx0 + int(bar_cols_rel[-1])
            matched = matched[:, bar_cols_rel[0] : bar_cols_rel[-1] + 1]
            nearest = nearest[:, bar_cols_rel[0] : bar_cols_rel[-1] + 1]
        else:
            bx0_ref, bx1_ref = bx0, bx1

        # Vectorized majority matched class per row across non-black columns
        min_row_hits = max(5, int(0.15 * matched.shape[1]))
        counts_per_row = np.stack([(matched & (nearest == k)).sum(axis=1) for k in range(n_cls)], axis=1)
        row_cls = np.where(matched.sum(axis=1) >= min_row_hits, np.argmax(counts_per_row, axis=1), -1)

        # Compute clean contiguous extent of each class k = 0..n_cls-1 (bottom to top)
        # Note: rows grow downward, so k=0 is near the bottom (large r_i) and k=n_cls-1 is near the top (small r_i).
        raw_counts = [int((row_cls == k).sum()) for k in range(n_cls)]
        tot_matched = sum(raw_counts)
        if tot_matched == 0:
            raise RuntimeError(f"{label} row={row_idx} bar={bi}: no bar pixels matched palette")

        # Measure exact boundaries from top/bottom row spans of each class, normalizing minor dark border pixels
        # by proportional allocation (so sum of all classes = 1.000 exactly and thin segments aren't biased by borders!)
        fracs = [c / float(tot_matched) for c in raw_counts]
        cum_tops = []
        acc = 0.0
        for f in fracs:
            acc += f
            cum_tops.append(min(1.0, acc))

        val = fracs[cls_idx]
        cum_lo = cum_tops[cls_idx - 1] if cls_idx > 0 else 0.0
        cum_hi = cum_tops[cls_idx]

        # Reading uncertainty from border/quantization (~2.5 px + calibration misfit)
        base_half_px = 2.5 + max(ya.residual_px, 0.5)
        half_val = base_half_px / px_per_unit

        sd_val = None
        stem_col = cx
        if has_sd:
            sd_val, stem_col, _ = _measure_staggered_sd(
                img, bx0_ref, bx1_ref, r0, r1, ya, cum_tops, fracs, cls_idx, n_cls, stem_fracs=stem_fracs
            )

        # Visual overlay marker sits right inside this bar at the segment's top boundary (or midpoint for top class)
        if cls_idx == n_cls - 1:
            # Place cross on the lower boundary of the top class (1 - val) where it meets class n_cls-2
            # (or at top boundary if val == 1.0) so it's clearly visible inside the bar!
            px_row = float(ya.to_pixel(cum_lo if val > 0.005 else 1.0))
        else:
            px_row = float(ya.to_pixel(cum_hi))

        lo_v = max(0.0, val - half_val)
        hi_v = min(1.0, val + half_val)
        rec: dict = {
            "x": float(xv),
            "y": float(val),
            "lo": float(lo_v),
            "hi": float(hi_v),
            "px": [float(stem_col), float(px_row)],
            "px_lo": float(px_row - base_half_px),
            "px_hi": float(px_row + base_half_px),
            "note": f"{label} row={row_idx} bar={bi+1} cls={cls_idx}: cum=[{cum_lo:.3f},{cum_hi:.3f}], ya_res={ya.residual_px:.2f}px",
        }
        if has_sd and sd_val is not None:
            sd_err = (2.0 + max(ya.residual_px, 0.5)) / px_per_unit
            b_lo = max(0.0, val - sd_val)
            b_hi = min(1.0, val + sd_val)
            rec["bar_lo"] = float(b_lo)
            rec["bar_hi"] = float(b_hi)
            rec["bar_lo_outer"] = float(max(0.0, b_lo - sd_err))
            rec["bar_hi_outer"] = float(min(1.0, b_hi + sd_err))
            # Draw the SD span right on the boundary where the black SD stem sits!
            if cls_idx == n_cls - 1:
                rec["px"] = [float(stem_col), float(ya.to_pixel(1.0 - sd_val))] if sd_val > 0.003 else [float(stem_col), float(px_row)]
                rec["px_lo"] = float(row_1)
                rec["px_hi"] = float(ya.to_pixel(max(0.0, 1.0 - sd_val)))
            else:
                rec["px_lo"] = float(ya.to_pixel(min(1.0, cum_hi + sd_val)))
                rec["px_hi"] = float(ya.to_pixel(max(0.0, cum_hi - sd_val)))
            rec["note"] += f"; SD={sd_val:.4f}"
        out.append(rec)
    return out
