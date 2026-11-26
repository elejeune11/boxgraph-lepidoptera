"""Shared calibrated reader for log-x marker+error-bar plots (Main Figs. 10B and 12B)."""
from __future__ import annotations

import math
import numpy as np
from scipy import ndimage


def _clusters(idx: np.ndarray, gap: int = 8) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap)[0] + 1
    return list(np.split(idx, splits))


def read_logx_red_markers(
    image: np.ndarray,
    at: list[float],
    y_ticks_vals: list[float],
    x_major_vals: list[float] = (0.1, 1.0, 10.0),
) -> list[dict]:
    """Calibrate log10-x / linear-y box and extract red circular markers + vertical error bars."""
    H, W = image.shape[:2]
    dark = (image[:, :, 0] < 80) & (image[:, :, 1] < 80) & (image[:, :, 2] < 80)

    # 1. Detect outer plot box borders
    h_rows = np.where(dark.sum(axis=1) > 0.45 * W)[0]
    hr = [float(c.mean()) for c in _clusters(h_rows, gap=10)]
    r0, r1 = hr[0], hr[-1]
    v_cols = np.where(dark[int(r0):int(r1), :].sum(axis=0) > 0.85 * (r1 - r0))[0]
    vc = [float(c.mean()) for c in _clusters(v_cols, gap=10)]
    c0, c1 = vc[0], vc[-1]

    # 2. Detect inward y ticks along left border and major x ticks along bottom border
    yt_mask = dark[int(r0 + 15):int(r1 - 15), int(c0 + 12):int(c0 + 35)].mean(axis=1) > 0.5
    yt_idx = np.where(yt_mask)[0] + int(r0 + 15)
    y_int = [float(c.mean()) for c in _clusters(yt_idx, gap=10)]
    y_tick_rows = np.array([r0] + y_int + [r1], dtype=float)
    y_vals_arr = np.array(y_ticks_vals, dtype=float)
    if len(y_tick_rows) != len(y_vals_arr):
        raise ValueError(f"Detected {len(y_tick_rows)} y-ticks, expected {len(y_vals_arr)}")
    yb, ya = np.polyfit(y_tick_rows, y_vals_arr, 1)
    y_res_px = float(np.max(np.abs(y_tick_rows - (y_vals_arr - ya) / yb)))

    # Major x-ticks extend farther up than minor ticks (25..45 px above r1)
    xt_mask = dark[int(r1 - 45):int(r1 - 25), int(c0 + 15):int(c1 - 15)].mean(axis=0) > 0.5
    xt_idx = np.where(xt_mask)[0] + int(c0 + 15)
    x_int = [float(c.mean()) for c in _clusters(xt_idx, gap=10)]
    x_tick_cols = np.array([c0] + x_int + [c1], dtype=float)
    x_log_arr = np.log10(np.array(x_major_vals, dtype=float))
    if len(x_tick_cols) != len(x_log_arr):
        raise ValueError(f"Detected {len(x_tick_cols)} x-ticks, expected {len(x_log_arr)}")
    xb, xa = np.polyfit(x_tick_cols, x_log_arr, 1)
    x_res_px = float(np.max(np.abs(x_tick_cols - (x_log_arr - xa) / xb)))

    # 3. Red mask for markers, error bars, and caps
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    red = (R > 150) & (R - G > 90) & (R - B > 90)

    # Structuring element (disk radius 12 px) isolates the ~45 px diameter circular marker interior,
    # completely suppressing ~6 px connecting lines and error-bar stems/caps.
    yy, xx = np.ogrid[-12:13, -12:13]
    disk = (xx * xx + yy * yy) <= 12 * 12

    xs = list(at) if at else [0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
    px_per_y = abs(1.0 / yb)
    out = []

    for xv in xs:
        cx_exp = (math.log10(float(xv)) - xa) / xb
        ci = int(round(cx_exp))
        strip_c0 = max(int(c0 + 4), ci - 45)
        strip_c1 = min(int(c1 - 4), ci + 46)
        strip_r0 = int(r0 + 6)
        strip_r1 = int(r1 - 6)

        sub_red = red[strip_r0:strip_r1, strip_c0:strip_c1]
        opened = ndimage.binary_opening(sub_red, structure=disk)
        lbl, n_comp = ndimage.label(opened)
        if n_comp == 0:
            out.append({"x": float(xv), "y": None, "note": "marker not found"})
            continue

        # Choose the connected component closest to cx_exp
        best_c, best_r, best_dist, best_r_half = None, None, 1e9, 20.0
        for cid in range(1, n_comp + 1):
            rr, cc = np.nonzero(lbl == cid)
            mc = float(cc.mean()) + strip_c0
            mr = float(rr.mean()) + strip_r0
            if abs(mc - cx_exp) < best_dist:
                best_dist = abs(mc - cx_exp)
                best_c, best_r = mc, mr
                best_r_half = 0.5 * (float(rr.max() - rr.min()) + 1.0)
        # Re-refine marker center using the full circular red pixels within radius+4 of (best_c, best_r)
        mr_i, mc_i = int(round(best_r)), int(round(best_c))
        rad_win = int(round(best_r_half + 10))
        wr0, wr1 = max(strip_r0, mr_i - rad_win), min(strip_r1, mr_i + rad_win + 1)
        wc0, wc1 = max(strip_c0, mc_i - rad_win), min(strip_c1, mc_i + rad_win + 1)
        # Black outline ring around the red circle gives exact outer radius and subpixel center
        ring_dark = dark[wr0:wr1, wc0:wc1]
        ring_red = red[wr0:wr1, wc0:wc1]
        # Dilate opened marker slightly to find exact filled circle pixels (excludes connecting line!)
        yy2, xx2 = np.ogrid[-(mr_i - wr0):(wr1 - mr_i), -(mc_i - wc0):(wc1 - mc_i)]
        near_disk = (xx2 * xx2 + yy2 * yy2) <= (best_r_half + 4.0) ** 2
        m_rr, m_cc = np.nonzero(ring_red & near_disk)
        cy = float(m_rr.mean()) + wr0
        cx = float(m_cc.mean()) + wc0
        # Marker outer radius (including black outline)
        r_outer = best_r_half + 5.0

        # Detect vertical error bar caps above and below the marker around column cx
        ci_m = int(round(cx))
        cap_c0 = max(int(c0 + 4), ci_m - 26)
        cap_c1 = min(int(c1 - 4), ci_m + 27)
        cap_strip = red[strip_r0:strip_r1, cap_c0:cap_c1]
        # Horizontal caps have >= 24 red pixels across the 53-px strip, and sit outside the marker!
        row_counts = cap_strip.sum(axis=1)
        # Also verify vertical red stem connects cap to marker at column ci_m +- 4
        stem_cols = red[strip_r0:strip_r1, max(0, ci_m - 4):min(W, ci_m + 5)].any(axis=1)

        # Find contiguous vertical stem extent starting from the marker center
        rel_cy = int(round(cy)) - strip_r0
        # Walk up while stem_cols or red/dark marker is present (allow small 4px gap across black outline)
        r_up = rel_cy
        gap_cnt = 0
        for rr_i in range(rel_cy, -1, -1):
            in_marker = abs(rr_i - rel_cy) <= r_outer + 3
            if stem_cols[rr_i] or in_marker:
                r_up = rr_i
                gap_cnt = 0
            else:
                gap_cnt += 1
                if gap_cnt > 5:
                    break
        r_dn = rel_cy
        gap_cnt = 0
        for rr_i in range(rel_cy, len(stem_cols)):
            in_marker = abs(rr_i - rel_cy) <= r_outer + 3
            if stem_cols[rr_i] or in_marker:
                r_dn = rr_i
                gap_cnt = 0
            else:
                gap_cnt += 1
                if gap_cnt > 5:
                    break

        top_px = float(r_up + strip_r0)
        bot_px = float(r_dn + strip_r0)
        hidden = (cy - top_px <= r_outer + 5.0) and (bot_px - cy <= r_outer + 5.0)

        if not hidden:
            # Refine cap center rows within 15 px of top_px and bot_px
            t_lo = max(0, r_up - 2)
            t_hi = min(len(row_counts), r_up + 18)
            t_rows = np.where(row_counts[t_lo:t_hi] >= 20)[0] + t_lo + strip_r0
            if len(t_rows):
                top_px = float(t_rows.mean())
            b_lo = max(0, r_dn - 18)
            b_hi = min(len(row_counts), r_dn + 3)
            b_rows = np.where(row_counts[b_lo:b_hi] >= 20)[0] + b_lo + strip_r0
            if len(b_rows):
                bot_px = float(b_rows.mean())
            note = f"cal_res=({x_res_px:.2f},{y_res_px:.2f})px"
        else:
            top_px = cy - r_outer
            bot_px = cy + r_outer
            note = f"error bar hidden inside marker (<= marker radius); cal_res=({x_res_px:.2f},{y_res_px:.2f})px"

        y_val = float(ya + yb * cy)
        y_u = abs(yb) * max(2.5, y_res_px + 1.5)
        bar_hi = float(ya + yb * top_px)
        bar_lo = float(ya + yb * bot_px)
        bar_lo, bar_hi = min(bar_lo, bar_hi), max(bar_lo, bar_hi)
        bar_u = abs(yb) * max(3.0, y_res_px + 2.0)

        out.append({
            "x": float(xv),
            "y": round(y_val, 5),
            "lo": round(y_val - y_u, 5),
            "hi": round(y_val + y_u, 5),
            "bar_lo": round(bar_lo, 5),
            "bar_hi": round(bar_hi, 5),
            "bar_lo_outer": round(bar_lo - bar_u, 5),
            "bar_hi_outer": round(bar_hi + bar_u, 5),
            "px": [round(cx, 2), round(cy, 2)],
            "px_lo": round(top_px, 2),
            "px_hi": round(bot_px, 2),
            "note": note,
        })
    return out
