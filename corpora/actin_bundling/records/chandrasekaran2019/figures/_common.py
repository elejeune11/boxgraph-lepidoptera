"""Shared fast, calibrated figure-reading helpers for chandrasekaran2019 (900-1200 dpi).

All helpers operate purely on the passed RGB image array (no PDF rendering or file I/O)
and support MATLAB-style plot frames with inward ticks on 2-spine (L-shaped) and 4-spine
rectangular axes.
"""

from __future__ import annotations

import numpy as np
from boxgraph import digitize as D


def _longest_true_run_1d(arr: np.ndarray, gap: int = 3) -> tuple[int, int, int]:
    """Longest True run in a 1D boolean array, bridging gaps <= gap."""
    idx = np.flatnonzero(arr)
    if idx.size == 0:
        return 0, 0, 0
    splits = np.flatnonzero(np.diff(idx) > (gap + 1))
    starts = np.concatenate(([0], splits + 1))
    ends = np.concatenate((splits, [idx.size - 1]))
    lengths = idx[ends] - idx[starts] + 1
    k = int(np.argmax(lengths))
    return int(lengths[k]), int(idx[starts[k]]), int(idx[ends[k]] + 1)


def find_spines_1d(mask_2d: np.ndarray, min_len: float, gap: int = 4, merge_dist: int = 5) -> list[dict]:
    """Find horizontal spine bands in a 2D boolean mask (rows with a long dark run >= min_len)."""
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


def find_ticks_in_strip(strip: np.ndarray, start_pos: int, end_pos: int,
                        min_len: int = 6, max_len: int = 45, max_width: int = 24) -> list[float]:
    """Find tick center positions along a boolean strip of shape (depth, L),
    where row 0 is immediately adjacent to the spine."""
    s = max(0, int(start_pos))
    e = min(strip.shape[1], int(end_pos))
    sub = strip[:, s:e]
    depth, width = sub.shape
    if depth == 0 or width == 0:
        return []
    # A tick starts at depth 0..2 and ends cleanly within max_len (so curves passing through are rejected)
    starts_ok = sub[0] | (sub[1] if depth > 1 else False) | (sub[2] if depth > 2 else False)
    lens = np.zeros(width, dtype=int)
    for k in np.flatnonzero(starts_ok):
        col = sub[:, k]
        j = 0
        while j < min(3, depth) and not col[j]:
            j += 1
        n = 0
        while j + n < depth and col[j + n]:
            n += 1
        lens[k] = n
    pos = np.flatnonzero((lens >= min_len) & (lens <= max_len))
    if pos.size == 0:
        return []
    splits = np.flatnonzero(np.diff(pos) > 2)
    g_starts = np.concatenate(([0], splits + 1))
    g_ends = np.concatenate((splits, [pos.size - 1]))
    ticks = []
    for gs, ge in zip(g_starts, g_ends):
        p0, p1 = int(pos[gs]), int(pos[ge])
        if (p1 - p0 + 1) <= max_width:
            ticks.append(s + 0.5 * (p0 + p1) + 0.5)
    return ticks


def find_matlab_box(img: np.ndarray, dark_thresh: float = 140.0,
                    min_h_frac: float = 0.45, min_v_frac: float = 0.45,
                    has_top: bool = True, has_right: bool = True,
                    min_tick_len: int = 6, max_tick_len: int = 45,
                    max_tick_w: int = 24, margin_px: int = 6,
                    neutral_only: bool = True) -> dict:
    """Locate a MATLAB rectangular (4-spine) or L-shaped (2-spine) plot frame and inward ticks."""
    r = img[:, :, 0]
    g = img[:, :, 1]
    b = img[:, :, 2]
    gray3 = r.astype(np.uint16) + g + b
    dark = gray3 < int(round(3.0 * dark_thresh))
    if neutral_only:
        cmax = np.maximum(np.maximum(r, g), b)
        cmin = np.minimum(np.minimum(r, g), b)
        dark &= (cmax - cmin) < 45
    H, W = dark.shape

    h_spines = find_spines_1d(dark, min_len=min_h_frac * W)
    v_spines = find_spines_1d(dark.T, min_len=min_v_frac * H)
    if not h_spines or not v_spines:
        raise ValueError(f"Could not locate plot spines (h={len(h_spines)}, v={len(v_spines)})")

    if has_top:
        # Filter to the pair of longest horizontal spines separated by > 25% of H
        top_sp = h_spines[0]
        bot_sp = h_spines[-1]
        if bot_sp["lo"] - top_sp["hi"] < 0.25 * H:
            raise ValueError(f"Top and bottom spines too close: {top_sp} vs {bot_sp}")
    else:
        bot_sp = max(h_spines, key=lambda s: s["length"])
        top_sp = None

    if has_right:
        left_sp = v_spines[0]
        right_sp = v_spines[-1]
        if right_sp["lo"] - left_sp["hi"] < 0.25 * W:
            raise ValueError(f"Left and right spines too close: {left_sp} vs {right_sp}")
    else:
        left_sp = max(v_spines, key=lambda s: s["length"])
        right_sp = None

    c0 = left_sp["hi"] + 1
    c1 = right_sp["lo"] if right_sp else bot_sp["end"]
    r0 = (top_sp["hi"] + 1) if top_sp else left_sp["start"]
    r1 = bot_sp["lo"]

    depth = max_tick_len + 8
    # Bottom inward ticks (pointing UP from bottom spine)
    x_bot_strip = dark[max(0, bot_sp["lo"] - depth):bot_sp["lo"], :][::-1, :]
    x_ticks_bot = find_ticks_in_strip(
        x_bot_strip, c0 + margin_px, c1 - margin_px,
        min_len=min_tick_len, max_len=max_tick_len, max_width=max_tick_w,
    )
    # Top inward ticks (pointing DOWN from top spine)
    x_ticks_top: list[float] = []
    if top_sp is not None:
        x_top_strip = dark[top_sp["hi"] + 1:min(H, top_sp["hi"] + 1 + depth), :]
        x_ticks_top = find_ticks_in_strip(
            x_top_strip, c0 + margin_px, c1 - margin_px,
            min_len=min_tick_len, max_len=max_tick_len, max_width=max_tick_w,
        )
    # Left inward ticks (pointing RIGHT from left spine)
    y_left_strip = dark[:, left_sp["hi"] + 1:min(W, left_sp["hi"] + 1 + depth)].T
    y_ticks_left = find_ticks_in_strip(
        y_left_strip, r0 + margin_px, r1 - margin_px,
        min_len=min_tick_len, max_len=max_tick_len, max_width=max_tick_w,
    )
    # Right inward ticks (pointing LEFT from right spine)
    y_ticks_right: list[float] = []
    if right_sp is not None:
        y_right_strip = dark[:, max(0, right_sp["lo"] - depth):right_sp["lo"]][:, ::-1].T
        y_ticks_right = find_ticks_in_strip(
            y_right_strip, r0 + margin_px, r1 - margin_px,
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
        "x_ticks_bot": x_ticks_bot,
        "x_ticks_top": x_ticks_top,
        "y_ticks_left": y_ticks_left,
        "y_ticks_right": y_ticks_right,
    }


def merge_tick_set(corner_lo: float, inner_ticks: list[float], corner_hi: float,
                   expected_n: int, tol_px: float = 14.0) -> list[float]:
    """Assemble an evenly-spaced tick set spanning [corner_lo, corner_hi] using detected inner ticks."""
    expected = np.linspace(corner_lo, corner_hi, expected_n)
    pts = [corner_lo]
    for target in expected[1:-1]:
        near = [t for t in inner_ticks if abs(t - target) <= tol_px]
        if near:
            pts.append(min(near, key=lambda t: abs(t - target)))
        else:
            raise ValueError(
                f"Missing inner tick near {target:.1f} px (inner={[round(t, 1) for t in inner_ticks]})"
            )
    pts.append(corner_hi)
    return pts


def fit_axis(ticks: list[float], values: list[float], scale: str = "linear",
             max_residual_px: float = 3.0) -> D.Axis:
    if len(ticks) != len(values):
        raise ValueError(
            f"Tick count mismatch: got {len(ticks)} ({[round(t, 1) for t in ticks]}), expected {len(values)} ({values})"
        )
    ax = D.Axis(scale, list(zip(ticks, values)))
    if ax.residual_px > max_residual_px:
        raise ValueError(f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px")
    return ax


def make_reading(x: float, row_px: float, y_axis: D.Axis, col_px: float, half_px: float,
                 x_axis: D.Axis | None = None, slope_px: float = 0.0, note: str = "",
                 bar_top_px: float | None = None, bar_bot_px: float | None = None,
                 bar_err_px: float = 2.0, clamp_y: tuple[float, float] | None = None) -> dict:
    """Build a BoxGraph mode:code reading dict with calibrated 95% uncertainty interval."""
    y_misfit = max(y_axis.residual_px, 0.75)
    x_misfit = max(x_axis.residual_px, 0.75) if x_axis is not None else 0.0
    h_px = float(half_px) + 0.75 + y_misfit + abs(slope_px) * x_misfit
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
        d["bar_lo"] = float(min(bar_lo, y_val))
        d["bar_hi"] = float(max(bar_hi, y_val))
        d["bar_lo_outer"] = float(min(bar_lo_outer, d["bar_lo"]))
        d["bar_hi_outer"] = float(max(bar_hi_outer, d["bar_hi"]))
    return d
