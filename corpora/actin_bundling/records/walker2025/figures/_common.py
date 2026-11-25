"""Shared calibrated figure-reading utilities for walker2025 (900-1200 dpi).

Every function calibrates plot frames and tick marks directly from pixel data in the
cropped figure image, computes honest 95% intervals including axis-calibration misfit,
line/edge width, and pixel quantization, and returns BoxGraph code-mode reading dicts.
"""

from __future__ import annotations

import math
import numpy as np
from boxgraph import digitize as D


def longest_true_run_1d(arr: np.ndarray, gap: int = 3) -> tuple[int, int, int]:
    """Longest True run in a 1D boolean array, bridging gaps <= gap.
    Returns (length, start_inclusive, end_exclusive)."""
    idx = np.flatnonzero(arr)
    if idx.size == 0:
        return 0, 0, 0
    splits = np.flatnonzero(np.diff(idx) > (gap + 1))
    starts = np.concatenate(([0], splits + 1))
    ends = np.concatenate((splits, [idx.size - 1]))
    lengths = idx[ends] - idx[starts] + 1
    k = int(np.argmax(lengths))
    return int(lengths[k]), int(idx[starts[k]]), int(idx[ends[k]] + 1)


def clusters_1d(idx: np.ndarray, gap: int = 2) -> list[tuple[float, int, int]]:
    """Group sorted 1D integer indices separated by <= gap.
    Returns list of (center, start_inclusive, end_inclusive) with continuous center."""
    if len(idx) == 0:
        return []
    arr = np.sort(np.asarray(idx, dtype=int))
    splits = np.flatnonzero(np.diff(arr) > gap)
    starts = np.concatenate(([0], splits + 1))
    ends = np.concatenate((splits, [len(arr) - 1]))
    out = []
    for s, e in zip(starts, ends):
        p0, p1 = int(arr[s]), int(arr[e])
        out.append((0.5 * (p0 + p1) + 0.5, p0, p1))
    return out


def group_spines(mask_2d: np.ndarray, min_len: float, gap: int = 4, merge_dist: int = 4) -> list[dict]:
    """Find horizontal spine bands in a 2D boolean mask (rows containing a run >= min_len)."""
    row_sums = mask_2d.sum(axis=1)
    cands = np.flatnonzero(row_sums >= 0.7 * min_len)
    hits = []
    for r in cands:
        length, s, e = longest_true_run_1d(mask_2d[r], gap=gap)
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


def find_ticks_on_strip(
    strip: np.ndarray,
    start_pos: int,
    end_pos: int,
    min_len: int = 4,
    max_len: int = 80,
    max_width: int = 24,
) -> list[float]:
    """Find tick center positions along a boolean strip of shape (depth, span),
    where row 0 is immediately adjacent to the axis spine."""
    s = max(0, int(start_pos))
    e = min(strip.shape[1], int(end_pos))
    sub = strip[:, s:e]
    depth, width = sub.shape
    if depth == 0 or width == 0:
        return []
    starts_ok = sub[0] | (sub[1] if depth > 1 else False) | (sub[2] if depth > 2 else False)
    lens = np.zeros(width, dtype=int)
    for k in np.flatnonzero(starts_ok):
        col = sub[:, k]
        j = 0
        while j < 3 and j < depth and not col[j]:
            j += 1
        n = 0
        while j + n < depth and col[j + n]:
            n += 1
        lens[k] = n
    pos = np.flatnonzero((lens >= min_len) & (lens <= max_len))
    if pos.size == 0:
        return []
    out = []
    for c, p0, p1 in clusters_1d(pos, gap=1):
        if (p1 - p0 + 1) <= max_width:
            out.append(s + c)
    return out


def match_regular_ticks(cands: list[float], n_expected: int, tol: float = 4.0) -> list[float]:
    """Select the best evenly-spaced subset of exactly n_expected tick positions from candidates.
    Prefers widest span with minimum residual."""
    pts = sorted(cands)
    if len(pts) < n_expected:
        raise ValueError(f"Only {len(pts)} tick candidates for {n_expected} expected: {pts}")
    if len(pts) == n_expected:
        return pts
    best = None
    best_score = None
    for i in range(len(pts)):
        for j in range(i + n_expected - 1, len(pts)):
            span = pts[j] - pts[i]
            if span < 30.0:
                continue
            step = span / (n_expected - 1)
            picked = []
            max_err = 0.0
            ok = True
            for k in range(n_expected):
                target = pts[i] + k * step
                near = min(pts, key=lambda x: abs(x - target))
                err = abs(near - target)
                if err > max(tol, 0.03 * step):
                    ok = False
                    break
                picked.append(near)
                max_err = max(max_err, err)
            if ok and len(set(picked)) == n_expected:
                # Prefer smaller residual, breaking near-ties by larger span
                score = (round(max_err, 2), -span)
                if best_score is None or score < best_score:
                    best_score = score
                    best = picked
    if best is None:
        raise ValueError(f"Could not find {n_expected} evenly spaced ticks among {pts}")
    return best


def find_l_axes(
    img: np.ndarray,
    dark_thresh: float = 120.0,
    min_h_frac: float = 0.35,
    min_v_frac: float = 0.35,
    min_tick_len: int = 4,
    max_tick_len: int = 60,
    max_x_tick_len: int | None = None,
    max_y_tick_len: int | None = None,
    max_tick_w: int = 32,
    margin_px: int = 6,
    boxed: bool = False,
) -> dict:
    """Locate an L-shaped or boxed plot frame (bottom x-axis, left y-axis, optional top/right)
    and outward ticks along bottom and left spines."""
    gray = img.astype(np.float64).mean(axis=2)
    dark = gray < dark_thresh
    H, W = dark.shape
    h_spines = group_spines(dark, min_len=min_h_frac * W)
    v_spines = group_spines(dark.T, min_len=min_v_frac * H)
    if not h_spines or not v_spines:
        raise ValueError(f"Could not find axes: h={len(h_spines)}, v={len(v_spines)}")
    bot = h_spines[-1]
    left = v_spines[0]
    top = (
        h_spines[0]
        if boxed and len(h_spines) >= 2 and h_spines[0]["hi"] < bot["lo"] - 20
        else None
    )
    right = (
        v_spines[-1]
        if boxed and len(v_spines) >= 2 and v_spines[-1]["lo"] > left["hi"] + 20
        else None
    )

    c0 = left["hi"] + 1
    c1 = right["lo"] if right else bot["end"]
    r0 = (top["hi"] + 1) if top else left["start"]
    r1 = bot["lo"]

    mx = max_x_tick_len if max_x_tick_len is not None else max_tick_len
    my = max_y_tick_len if max_y_tick_len is not None else max_tick_len
    x_strip = dark[bot["hi"] + 1 : min(H, bot["hi"] + mx + 6), :]
    x_ticks = find_ticks_on_strip(
        x_strip,
        left["lo"] - margin_px,
        (right["hi"] if right else bot["end"]) + margin_px + 1,
        min_len=min_tick_len,
        max_len=mx,
        max_width=max_tick_w,
    )
    y_strip = dark[:, max(0, left["lo"] - my - 6) : left["lo"]][:, ::-1].T
    y_ticks = find_ticks_on_strip(
        y_strip,
        (top["lo"] if top else left["start"]) - margin_px,
        bot["hi"] + margin_px + 1,
        min_len=min_tick_len,
        max_len=my,
        max_width=max_tick_w,
    )
    return {
        "box": (int(c0), int(r0), int(c1), int(r1)),
        "bot": bot,
        "left": left,
        "top": top,
        "right": right,
        "x_ticks": x_ticks,
        "y_ticks": y_ticks,
    }


def make_axis(scale: str, marks: list[tuple[float, float]], max_res_px: float = 4.0) -> D.Axis:
    ax = D.Axis(scale, marks)
    if ax.residual_px > max_res_px:
        raise ValueError(f"Axis calibration residual {ax.residual_px:.2f} px > {max_res_px} px; marks={marks}")
    return ax


def reading_point(
    x_val: float,
    col: float,
    row: float,
    y_axis: D.Axis,
    half_px: float,
    note: str = "",
    bar_hi_row: float | None = None,
    bar_lo_row: float | None = None,
    bar_err_px: float | None = None,
    clamp_min: float | None = None,
    clamp_max: float | None = None,
) -> dict:
    """Build a BoxGraph code-mode reading dict with honest calibration+width uncertainty."""
    h = max(half_px, 0.5) + max(y_axis.residual_px, 0.5)
    y = y_axis.to_value(row)
    a, b = y_axis.to_value(row - h), y_axis.to_value(row + h)
    lo, hi = min(a, b), max(a, b)
    if clamp_min is not None:
        if y < clamp_min and hi >= clamp_min:
            y = clamp_min
        lo = max(clamp_min, lo)
        hi = max(y, hi)
    if clamp_max is not None:
        if y > clamp_max and lo <= clamp_max:
            y = clamp_max
        hi = min(clamp_max, hi)
        lo = min(y, lo)
    rec: dict = {
        "x": float(x_val),
        "y": float(y),
        "lo": float(lo),
        "hi": float(hi),
        "px": [float(col), float(row)],
        "px_lo": float(row - h),
        "px_hi": float(row + h),
        "note": note,
    }
    if bar_hi_row is not None and bar_lo_row is not None:
        be = (bar_err_px if bar_err_px is not None else 1.5) + max(y_axis.residual_px, 0.5)
        # Note: smaller pixel row = larger y value
        r_top = min(bar_hi_row, bar_lo_row)
        r_bot = max(bar_hi_row, bar_lo_row)
        v_hi = y_axis.to_value(r_top)
        v_lo = y_axis.to_value(r_bot)
        v_hi_out = y_axis.to_value(r_top - be)
        v_lo_out = y_axis.to_value(r_bot + be)
        b_lo, b_hi = min(v_lo, v_hi), max(v_lo, v_hi)
        b_lo_o, b_hi_o = min(v_lo_out, v_hi_out), max(v_lo_out, v_hi_out)
        if clamp_min is not None:
            b_lo = max(clamp_min, b_lo)
            b_lo_o = max(clamp_min, b_lo_o)
        if clamp_max is not None:
            b_hi = min(clamp_max, b_hi)
            b_hi_o = min(clamp_max, b_hi_o)
        b_lo = min(b_lo, y)
        b_hi = max(b_hi, y)
        b_lo_o = min(b_lo_o, b_lo)
        b_hi_o = max(b_hi_o, b_hi)
        rec["bar_lo"] = float(b_lo)
        rec["bar_hi"] = float(b_hi)
        rec["bar_lo_outer"] = float(b_lo_o)
        rec["bar_hi_outer"] = float(b_hi_o)
        # Draw the error bar span on the overlay
        rec["px_lo"] = float(r_top)
        rec["px_hi"] = float(r_bot)
    return rec
