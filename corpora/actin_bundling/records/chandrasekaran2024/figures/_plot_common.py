"""Shared fast, calibrated figure-reading helpers for chandrasekaran2024 (900-1200 dpi)."""

from __future__ import annotations

import math
import numpy as np
from boxgraph import digitize as D


def _longest_true_run_1d(arr: np.ndarray, gap: int = 3) -> tuple[int, int, int]:
    """Fast vectorized longest True run in a 1D boolean array, bridging gaps <= gap."""
    idx = np.flatnonzero(arr)
    if idx.size == 0:
        return 0, 0, 0
    # Split where gap between consecutive True indices exceeds gap + 1
    splits = np.flatnonzero(np.diff(idx) > (gap + 1))
    starts = np.concatenate(([0], splits + 1))
    ends = np.concatenate((splits, [idx.size - 1]))
    lengths = idx[ends] - idx[starts] + 1
    k = int(np.argmax(lengths))
    return int(lengths[k]), int(idx[starts[k]]), int(idx[ends[k]] + 1)


def _group_spine_indices(mask_2d: np.ndarray, min_len: float, gap: int = 4, merge_dist: int = 4):
    """Find horizontal spine bands in a 2D boolean mask (rows with a long run >= min_len)."""
    # Fast prefilter: row must have at least 0.7 * min_len True pixels
    row_sums = mask_2d.sum(axis=1)
    cands = np.flatnonzero(row_sums >= 0.7 * min_len)
    hits = []
    for r in cands:
        length, s, e = _longest_true_run_1d(mask_2d[r], gap=gap)
        if length >= min_len:
            hits.append((int(r), length, s, e))
    groups = []
    for h in hits:
        if not groups or h[0] > groups[-1][-1][0] + merge_dist:
            groups.append([h])
        else:
            groups[-1].append(h)
    out = []
    for g in groups:
        r0, r1 = g[0][0], g[-1][0]
        best = max(g, key=lambda x: x[1])
        out.append({
            "lo": r0,
            "hi": r1,
            "center": 0.5 * (r0 + r1) + 0.5,
            "start": best[2],
            "end": best[3],
            "length": best[1],
        })
    return out


def _find_ticks_on_strip(strip: np.ndarray, start_pos: int, end_pos: int,
                         min_len: int = 5, max_len: int = 90, max_width: int = 28) -> list[float]:
    """Find tick center coordinates along a boolean strip of shape (depth, W_or_H)
    where index 0 along axis 0 is immediately adjacent to the spine."""
    s = max(0, int(start_pos))
    e = min(strip.shape[1], int(end_pos))
    sub = strip[:, s:e]
    depth, width = sub.shape
    if depth == 0 or width == 0:
        return []
    # Run length starting within the first 2 pixels from the spine
    # Vectorized: find first True in rows 0..2, then first False after it
    starts_ok = sub[0] | (sub[1] if depth > 1 else False)
    lens = np.zeros(width, dtype=int)
    for k in np.flatnonzero(starts_ok):
        col = sub[:, k]
        j = 0
        while j < 2 and j < depth and not col[j]:
            j += 1
        n = 0
        while j + n < depth and col[j + n]:
            n += 1
        lens[k] = n
    pos = np.flatnonzero((lens >= min_len) & (lens <= max_len))
    if pos.size == 0:
        return []
    splits = np.flatnonzero(np.diff(pos) > 1)
    g_starts = np.concatenate(([0], splits + 1))
    g_ends = np.concatenate((splits, [pos.size - 1]))
    ticks = []
    for gs, ge in zip(g_starts, g_ends):
        p0, p1 = int(pos[gs]), int(pos[ge])
        if (p1 - p0 + 1) <= max_width:
            ticks.append(s + 0.5 * (p0 + p1) + 0.5)
    return ticks


def find_plot_box(img: np.ndarray, dark_thresh: float = 110.0, min_frac: float = 0.42,
                  max_tick_w: int = 28, min_tick_len: int = 5, max_tick_len: int = 90,
                  margin_px: int = 2) -> dict:
    """Locate the main rectangular plot frame (left, top, right, bottom spines) and
    outward tick marks along the bottom x-axis and left y-axis."""
    gray = img.astype(np.float64).mean(axis=2)
    dark = gray < dark_thresh
    H, W = dark.shape
    h_spines = _group_spine_indices(dark, min_len=min_frac * W)
    v_spines = _group_spine_indices(dark.T, min_len=min_frac * H)
    if not h_spines or not v_spines:
        raise ValueError(f"Could not find plot spines (h={len(h_spines)}, v={len(v_spines)})")
    # Pick spines that form a coherent enclosing box:
    # Bottom x-axis is the lowest long horizontal spine; left y-axis is the leftmost long vertical spine.
    bot_sp = h_spines[-1]
    top_sp = h_spines[0] if len(h_spines) >= 2 and h_spines[0]["hi"] < bot_sp["lo"] - 20 else None
    left_sp = v_spines[0]
    right_sp = v_spines[-1] if len(v_spines) >= 2 and v_spines[-1]["lo"] > left_sp["hi"] + 20 else None

    c0 = left_sp["hi"] + 1
    c1 = right_sp["lo"] if right_sp else bot_sp["end"]
    r0 = (top_sp["hi"] + 1) if top_sp else left_sp["start"]
    r1 = bot_sp["lo"]

    x_strip = dark[bot_sp["hi"] + 1:min(H, bot_sp["hi"] + max_tick_len + 6), :]
    x_ticks = _find_ticks_on_strip(
        x_strip, left_sp["lo"] - margin_px, (right_sp["hi"] if right_sp else bot_sp["end"]) + margin_px + 1,
        min_len=min_tick_len, max_len=max_tick_len, max_width=max_tick_w,
    )
    y_strip = dark[:, max(0, left_sp["lo"] - max_tick_len - 6):left_sp["lo"]][:, ::-1].T
    y_ticks = _find_ticks_on_strip(
        y_strip, (top_sp["lo"] if top_sp else left_sp["start"]) - margin_px, bot_sp["hi"] + margin_px + 1,
        min_len=min_tick_len, max_len=max_tick_len, max_width=max_tick_w,
    )
    return {
        "box": (int(c0), int(r0), int(c1), int(r1)),
        "spines": {
            "left": left_sp["center"],
            "right": right_sp["center"] if right_sp else float(bot_sp["end"]),
            "top": top_sp["center"] if top_sp else float(left_sp["start"]),
            "bottom": bot_sp["center"],
            "left_band": (left_sp["lo"], left_sp["hi"]),
            "right_band": (right_sp["lo"], right_sp["hi"]) if right_sp else None,
            "top_band": (top_sp["lo"], top_sp["hi"]) if top_sp else None,
            "bottom_band": (bot_sp["lo"], bot_sp["hi"]),
        },
        "x_ticks": x_ticks,
        "y_ticks": y_ticks,
    }


def find_stacked_boxes(img: np.ndarray, n_rows: int, dark_thresh: float = 110.0,
                       min_h_frac: float = 0.40, min_v_frac: float = 0.08,
                       max_tick_w: int = 24, min_tick_len: int = 4, max_tick_len: int = 70) -> list[dict]:
    """Find `n_rows` vertically stacked subpanel boxes (e.g., Supp. Fig. 2C column, Fig. 3B/3D, Supp. Fig. 11D)."""
    gray = img.astype(np.float64).mean(axis=2)
    dark = gray < dark_thresh
    H, W = dark.shape
    h_spines = _group_spine_indices(dark, min_len=min_h_frac * W)
    # Pair consecutive horizontal spines into n_rows boxes (2*n_rows spines if separate, or n_rows+1 if shared)
    if len(h_spines) == 2 * n_rows:
        pairs = [(h_spines[2 * i], h_spines[2 * i + 1]) for i in range(n_rows)]
    elif len(h_spines) == n_rows + 1:
        pairs = [(h_spines[i], h_spines[i + 1]) for i in range(n_rows)]
    else:
        # Filter to spines with similar horizontal span as the widest spine
        max_len = max(s["length"] for s in h_spines)
        filt = [s for s in h_spines if s["length"] >= 0.85 * max_len]
        if len(filt) == 2 * n_rows:
            pairs = [(filt[2 * i], filt[2 * i + 1]) for i in range(n_rows)]
        elif len(filt) == n_rows + 1:
            pairs = [(filt[i], filt[i + 1]) for i in range(n_rows)]
        else:
            raise ValueError(f"Expected {2*n_rows} or {n_rows+1} h_spines, got {len(h_spines)} (filt={len(filt)})")
    out = []
    for top_sp, bot_sp in pairs:
        sub_dark = dark[top_sp["lo"]:bot_sp["hi"] + 1, :]
        sub_H = sub_dark.shape[0]
        v_sp = _group_spine_indices(sub_dark.T, min_len=0.75 * sub_H)
        left_sp, right_sp = v_sp[0], v_sp[-1]
        c0, c1 = left_sp["hi"] + 1, right_sp["lo"]
        r0, r1 = top_sp["hi"] + 1, bot_sp["lo"]
        x_strip = dark[bot_sp["hi"] + 1:min(H, bot_sp["hi"] + max_tick_len + 6), :]
        xt = _find_ticks_on_strip(x_strip, left_sp["lo"] - 2, right_sp["hi"] + 3,
                                  min_len=min_tick_len, max_len=max_tick_len, max_width=max_tick_w)
        y_strip = dark[:, max(0, left_sp["lo"] - max_tick_len - 6):left_sp["lo"]][:, ::-1].T
        yt = _find_ticks_on_strip(y_strip, top_sp["lo"] - 2, bot_sp["hi"] + 3,
                                  min_len=min_tick_len, max_len=max_tick_len, max_width=max_tick_w)
        out.append({
            "box": (int(c0), int(r0), int(c1), int(r1)),
            "spines": {
                "left": left_sp["center"],
                "right": right_sp["center"],
                "top": top_sp["center"],
                "bottom": bot_sp["center"],
            },
            "x_ticks": xt,
            "y_ticks": yt,
        })
    return out


def fit_axis(ticks: list[float], values: list[float], scale: str = "linear",
             max_residual_px: float = 2.0) -> D.Axis:
    if len(ticks) != len(values):
        raise ValueError(f"Tick count mismatch: found {len(ticks)} ticks ({[round(t,1) for t in ticks]}), expected {len(values)} values ({values})")
    ax = D.Axis(scale, list(zip(ticks, values)))
    if ax.residual_px > max_residual_px:
        raise ValueError(f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px")
    return ax


def make_reading(x: float, row_px: float, y_axis: D.Axis, col_px: float, half_px: float,
                 x_axis: D.Axis | None = None, slope_px: float = 0.0, note: str = "",
                 bar_top_px: float | None = None, bar_bot_px: float | None = None,
                 bar_err_px: float = 1.5, clamp_y: tuple[float, float] | None = None) -> dict:
    """Build a BoxGraph mode:code reading dict with honest 95% calibration+localization uncertainty."""
    y_misfit = max(y_axis.residual_px, 0.5)
    x_misfit = max(x_axis.residual_px, 0.5) if x_axis is not None else 0.0
    h_px = float(half_px) + 0.5 + y_misfit + abs(slope_px) * x_misfit
    y_val = y_axis.to_value(row_px)
    v1 = y_axis.to_value(row_px - h_px)
    v2 = y_axis.to_value(row_px + h_px)
    lo, hi = min(v1, v2), max(v1, v2)
    if clamp_y is not None:
        c_lo, c_hi = clamp_y
        y_val = min(max(y_val, c_lo), c_hi)
        lo = min(max(lo, c_lo), y_val)
        hi = max(min(hi, c_hi), y_val)
    d = {
        "x": float(x),
        "y": float(y_val),
        "lo": float(lo),
        "hi": float(hi),
        "px": [float(col_px), float(row_px)],
        "px_lo": float(row_px - h_px),
        "px_hi": float(row_px + h_px),
        "note": note,
    }
    if bar_top_px is not None and bar_bot_px is not None:
        e_px = float(bar_err_px) + y_misfit
        b1, b2 = y_axis.to_value(bar_top_px), y_axis.to_value(bar_bot_px)
        bar_lo, bar_hi = min(b1, b2), max(b1, b2)
        bo1, bo2 = y_axis.to_value(bar_top_px - e_px), y_axis.to_value(bar_bot_px + e_px)
        bar_lo_outer, bar_hi_outer = min(bo1, bo2), max(bo1, bo2)
        if clamp_y is not None:
            c_lo, c_hi = clamp_y
            bar_lo = max(bar_lo, c_lo)
            bar_hi = min(bar_hi, c_hi)
            bar_lo_outer = max(bar_lo_outer, c_lo)
            bar_hi_outer = min(bar_hi_outer, c_hi)
        d["bar_lo"] = float(bar_lo)
        d["bar_hi"] = float(bar_hi)
        d["bar_lo_outer"] = float(min(bar_lo_outer, bar_lo))
        d["bar_hi_outer"] = float(max(bar_hi_outer, bar_hi))
    return d


def rgb_dist(img: np.ndarray, rgb: tuple[float, float, float]) -> np.ndarray:
    """Euclidean RGB distance map."""
    target = np.asarray(rgb, dtype=np.float64)
    diff = img.astype(np.float64) - target
    return np.sqrt(np.sum(diff * diff, axis=2))


def sample_curve_columns(mask: np.ndarray, x_axis: D.Axis, y_axis: D.Axis, at: list[float],
                         search_half_w: int = 3, y_hint: dict[float, float] | None = None,
                         prefer: str = "median", min_px_per_col: int = 1,
                         clamp_y: tuple[float, float] | None = None,
                         note_prefix: str = "") -> list[dict]:
    """Read a single clean curve from a binary pixel mask at requested x values."""
    H, W = mask.shape
    out = []
    for xv in at:
        pc = x_axis.to_pixel(xv)
        col = int(round(pc))
        # Collect runs in small window around col
        best_col = None
        runs = []
        for dc in sorted(range(-search_half_w, search_half_w + 1), key=abs):
            c = col + dc
            if 0 <= c < W:
                rs = np.flatnonzero(mask[:, c])
                if rs.size >= min_px_per_col:
                    splits = np.flatnonzero(np.diff(rs) > 2)
                    st = np.concatenate(([0], splits + 1))
                    en = np.concatenate((splits, [rs.size - 1]))
                    runs = [(0.5 * (rs[a] + rs[b]) + 0.5, 0.5 * (rs[b] - rs[a] + 1)) for a, b in zip(st, en)]
                    best_col = c
                    break
        if not runs:
            out.append({"x": float(xv), "y": None, "note": f"no pixels found near x={xv:g}"})
            continue
        if len(runs) == 1:
            center, half = runs[0]
        elif y_hint is not None and xv in y_hint:
            target_r = y_axis.to_pixel(y_hint[xv])
            center, half = min(runs, key=lambda r: abs(r[0] - target_r))
        elif prefer == "top":
            center, half = min(runs, key=lambda r: r[0])
        elif prefer == "bottom":
            center, half = max(runs, key=lambda r: r[0])
        else:
            # Prefer widest run or nearest to local neighbor median
            neigh = []
            for dc in (-5, -3, 3, 5):
                c2 = best_col + dc
                if 0 <= c2 < W:
                    rs2 = np.flatnonzero(mask[:, c2])
                    if rs2.size:
                        neigh.append(float(np.median(rs2)) + 0.5)
            if neigh:
                ref = float(np.median(neigh))
                center, half = min(runs, key=lambda r: abs(r[0] - ref))
            else:
                center, half = max(runs, key=lambda r: r[1])
        # Estimate local slope for x-uncertainty propagation
        slope = 0.0
        left_r, right_r = None, None
        for d_off in (2, 3, 4):
            if best_col - d_off >= 0 and best_col + d_off < W:
                rl = np.flatnonzero(mask[:, best_col - d_off])
                rr = np.flatnonzero(mask[:, best_col + d_off])
                if rl.size and rr.size:
                    # Pick pixels close to center
                    rl_near = rl[np.abs(rl + 0.5 - center) < max(30.0, 4.0 * half)]
                    rr_near = rr[np.abs(rr + 0.5 - center) < max(30.0, 4.0 * half)]
                    if rl_near.size and rr_near.size:
                        left_r = float(np.median(rl_near)) + 0.5
                        right_r = float(np.median(rr_near)) + 0.5
                        slope = (right_r - left_r) / (2.0 * d_off)
                        break
        offset = abs(pc - (best_col + 0.5))
        eff_half = half + abs(slope) * offset
        note = f"{note_prefix} (rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)".strip()
        out.append(make_reading(xv, center, y_axis, pc, eff_half, x_axis=x_axis,
                                slope_px=slope, note=note, clamp_y=clamp_y))
    return out
