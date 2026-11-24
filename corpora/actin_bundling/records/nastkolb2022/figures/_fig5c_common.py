"""Shared calibration and tri-axis mean +/- SD curve reader for Fig. 5C (PDF p. 8, 1200 dpi)."""

from __future__ import annotations

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (
    LinearCal,
    detect_ticks_along_spine,
    select_major_ticks,
    split_runs,
)


def read_fig5c(image: np.ndarray, channel: str, at: list[float]) -> list[dict]:
    """Read Fig. 5C tri-axis mean +/- SD curves ('actin', 'vasp', 'lipids') at requested times (min).

    Detects true vertical error bars with horizontal T-caps, bridging crossings where another channel's
    curve passes over an error bar stem, and reads exact top/bottom SD caps and their midpoint.
    """
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]

    blk = (R < 85) & (G < 85) & (B < 85)
    red_ax = (R > 130) & (R - G > 70) & (R - B > 70)
    grn_ax = (G > 75) & (G - R > 25) & (G - B > 10)
    blu_ax = (B > 115) & (B - R > 60)

    # Spines
    hr_runs = split_runs(np.where(blk.sum(axis=1) > 0.35 * W)[0], gap=6)
    hr = (int(hr_runs[-1][0]), int(hr_runs[-1][-1]))
    vl_runs = split_runs(np.where(red_ax.sum(axis=0) > 0.40 * H)[0], gap=6)
    vl = (int(vl_runs[0][0]), int(vl_runs[0][-1]))
    vg_runs = split_runs(np.where(grn_ax.sum(axis=0) > 0.40 * H)[0], gap=6)
    vg = (int(vg_runs[-1][0]), int(vg_runs[-1][-1]))
    vb_runs = split_runs(np.where(blu_ax.sum(axis=0) > 0.40 * H)[0], gap=6)
    vb = (int(vb_runs[-1][0]), int(vb_runs[-1][-1]))

    # X axis: 0, 10, 20, 30, 40, 50, 60 min
    xt_raw = detect_ticks_along_spine(blk, hr, "x", "down", min_offset=2, max_offset=45, min_len=8, coord_min=vl[1], coord_max=vg[1] + 30)
    xt = select_major_ticks(xt_raw, expected=7, frac=0.70)
    x_cal = LinearCal(xt, [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], max_res_px=2.0)

    if channel == "actin":
        yt_raw = detect_ticks_along_spine(red_ax, vl, "y", "left", min_offset=2, max_offset=45, min_len=8, coord_min=50, coord_max=hr[0] - 5)
        yt = select_major_ticks(yt_raw, expected=8, frac=0.70)
        y_cal = LinearCal(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2, -0.4], max_res_px=2.0)
        cmask = (R > 135) & (R - G > 65) & (R - B > 65)
        # Mask out the "Actin" text label (top right, t > 52, y > 0.93)
        cmask[0 : int(y_cal.to_px(0.93)), int(x_cal.to_px(52.0)) : vg[0]] = False
        other_mask = grn_ax | blu_ax
        r_top = int(y_cal.to_px(1.06))
        r_bot = int(y_cal.to_px(-0.10))
        cap_half = 11
        stem_min_h = 42
    elif channel == "vasp":
        yt_raw = detect_ticks_along_spine(grn_ax, vg, "y", "right", min_offset=2, max_offset=45, min_len=8, coord_min=50, coord_max=hr[0] - 5)
        yt = select_major_ticks(yt_raw, expected=10, frac=0.70)
        y_cal = LinearCal(yt, [2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3], max_res_px=2.5)
        cmask = (G > 75) & (G - R > 22) & (G - B > 8)
        cmask[int(y_cal.to_px(1.94)) : int(y_cal.to_px(1.80)), int(x_cal.to_px(49.0)) : int(x_cal.to_px(59.5))] = False
        other_mask = red_ax | blu_ax
        r_top = int(y_cal.to_px(2.22))
        r_bot = hr[0] - 6
        cap_half = 11
        stem_min_h = 180
    elif channel == "lipids":
        yt_raw = detect_ticks_along_spine(blu_ax, vb, "y", "right", min_offset=2, max_offset=45, min_len=8, coord_min=50, coord_max=hr[0] - 5)
        yt = select_major_ticks(yt_raw, expected=8, frac=0.70)
        y_cal = LinearCal(yt, [1.25, 1.20, 1.15, 1.10, 1.05, 1.00, 0.95, 0.90], max_res_px=2.0)
        cmask = (B > 125) & (B - R > 65)
        cmask[int(y_cal.to_px(0.916)) : hr[0], int(x_cal.to_px(53.0)) : vg[0]] = False
        other_mask = red_ax | grn_ax
        r_top = int(y_cal.to_px(1.03))
        r_bot = hr[0] - 4
        cap_half = 10
        stem_min_h = 88
    else:
        raise ValueError(channel)

    c_left = vl[1] + 8
    c_right = vg[0] - 8
    sub = np.zeros((H, W), dtype=bool)
    sub[r_top:r_bot, c_left:c_right] = cmask[r_top:r_bot, c_left:c_right]

    # Detect error bars: require vertical stem (bridging gaps where another color overlaps) AND horizontal T-caps!
    bars = []  # (c_center, r_top_cap, r_bot_cap)
    # First find columns that have strong vertical presence in `sub`
    col_span = np.zeros(W, dtype=int)
    col_r0 = np.full(W, -1, dtype=int)
    col_r1 = np.full(W, -1, dtype=int)
    for c in range(c_left, c_right):
        own_rows = np.nonzero(sub[r_top:r_bot, c])[0] + r_top
        if len(own_rows) < 20:
            continue
        # Bridge gaps across other channel pixels (e.g. blue/red curves crossing green bar at t=0)
        bridged = sub[r_top:r_bot, c] | other_mask[r_top:r_bot, c]
        # Walk runs in `bridged` that both start and end on `own_rows`
        grs = split_runs(np.nonzero(bridged)[0] + r_top, gap=14)
        best_span = 0
        best_pair = None
        for gr in grs:
            own_in = gr[sub[gr, c]]
            if len(own_in) < 25:
                continue
            sp = int(own_in[-1] - own_in[0] + 1)
            if sp > best_span:
                best_span = sp
                best_pair = (int(own_in[0]), int(own_in[-1]))
        if best_pair is not None and best_span >= stem_min_h:
            r0, r1 = best_pair
            # Verify horizontal cap support at top AND bottom (at least 60% of 2*cap_half+1 width!)
            top_cap = sub[max(r_top, r0 - 3) : min(r_bot, r0 + 10), max(c_left, c - cap_half) : min(c_right, c + cap_half + 1)].any(axis=0).mean()
            bot_cap = sub[max(r_top, r1 - 10) : min(r_bot, r1 + 4), max(c_left, c - cap_half) : min(c_right, c + cap_half + 1)].any(axis=0).mean()
            if top_cap >= 0.65 and bot_cap >= 0.65:
                col_span[c] = best_span
                col_r0[c] = r0
                col_r1[c] = r1

    stem_groups = split_runs(np.where(col_span >= stem_min_h)[0], gap=4)
    for gr in stem_groups:
        if len(gr) <= 14:
            c_mid = float(gr.mean()) + 0.5
            r_hi = float(col_r0[gr].min()) + 0.5
            r_lo = float(col_r1[gr].max()) + 0.5
            bars.append((c_mid, r_hi, r_lo))
        else:
            # Wide merged cluster of adjacent bars (early dense times): pick local peaks
            for offset in range(0, len(gr), 10):
                chunk = gr[offset : offset + 10]
                bc = chunk[int(np.argmax(col_span[chunk]))]
                bars.append((float(bc) + 0.5, float(col_r0[bc]) + 0.5, float(col_r1[bc]) + 0.5))

    out = []
    px_per_min = abs(x_cal.to_px(10.0) - x_cal.to_px(0.0)) / 10.0

    for tv in at:
        pc = x_cal.to_px(tv)
        max_snap = (0.45 if tv <= 0.1 else 1.8) * px_per_min
        cand = [b for b in bars if abs(b[0] - pc) <= max_snap]
        if cand:
            best_b = min(cand, key=lambda b: abs(b[0] - pc))
            c_b, r_hi, r_lo = best_b
            ci = int(round(c_b))
            cap_band = sub[:, max(c_left, ci - cap_half) : min(c_right, ci + cap_half + 1)].mean(axis=1)
            top_rows = np.where(cap_band[max(r_top, int(r_hi) - 5) : min(r_bot, int(r_hi) + 12)] > 0.55)[0]
            if len(top_rows) > 0:
                r_hi = float(top_rows.mean()) + 0.5 + max(r_top, int(r_hi) - 5)
            bot_rows = np.where(cap_band[max(r_top, int(r_lo) - 12) : min(r_bot, int(r_lo) + 6)] > 0.55)[0]
            if len(bot_rows) > 0:
                r_lo = float(bot_rows.mean()) + 0.5 + max(r_top, int(r_lo) - 12)
            r_mid = 0.5 * (r_hi + r_lo)
            t_actual = x_cal.to_val(c_b)
            y_val = y_cal.to_val(r_mid)
            bar_hi = y_cal.to_val(r_hi)
            bar_lo = y_cal.to_val(r_lo)
            err_y = abs(y_cal.b) * (3.0 + y_cal.res_px)
            bar_pad = abs(y_cal.b) * (4.0 + y_cal.res_px)
            out.append({
                "x": round(t_actual, 2),
                "y": y_val,
                "lo": y_val - err_y,
                "hi": y_val + err_y,
                "bar_lo": min(bar_lo, bar_hi),
                "bar_hi": max(bar_lo, bar_hi),
                "bar_lo_outer": min(bar_lo, bar_hi) - bar_pad,
                "bar_hi_outer": max(bar_lo, bar_hi) + bar_pad,
                "px": [c_b, r_mid],
                "px_lo": r_hi,
                "px_hi": r_lo,
                "note": f"SD bar midpoint at t={t_actual:.2f}min (req {tv:g}), SD={0.5*abs(bar_hi-bar_lo):.4g}, rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
            })
        else:
            # Point without visible error bar (e.g. Actin at t=0.0 min where normalized value is 0.00)
            col = int(round(pc))
            rows = []
            for dc in range(-4, 5):
                c = col + dc
                if c_left <= c < c_right:
                    rs = np.nonzero(sub[r_top:r_bot, c])[0] + r_top
                    if len(rs) > 0:
                        rows.append(float(np.median(rs)) + 0.5)
            if not rows:
                out.append({"x": float(tv), "y": None, "found": False, "note": f"no pixels near t={tv:g}"})
                continue
            r_mid = float(np.median(rows))
            y_val = y_cal.to_val(r_mid)
            err_y = abs(y_cal.b) * (4.0 + y_cal.res_px)
            out.append({
                "x": float(tv),
                "y": y_val,
                "lo": y_val - err_y,
                "hi": y_val + err_y,
                "bar_lo": y_val - err_y,
                "bar_hi": y_val + err_y,
                "bar_lo_outer": y_val - 2 * err_y,
                "bar_hi_outer": y_val + 2 * err_y,
                "px": [pc, r_mid],
                "px_lo": r_mid - 4.0,
                "px_hi": r_mid + 4.0,
                "note": f"curve point without visible SD bar at t={tv:g}, rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
            })
    return out
