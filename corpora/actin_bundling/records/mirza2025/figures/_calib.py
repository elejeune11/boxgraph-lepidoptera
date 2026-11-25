"""Shared image-based axis calibration and curve/marker digitization helpers
for mirza2025. Every coordinate, tick position, curve point, marker centroid,
and error bar end is measured directly from the rendered image pixels.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from boxgraph import digitize as D


def groups(arr: np.ndarray, gap: int = 2) -> list[tuple[int, int]]:
    """Split a sorted 1D integer index array into contiguous (start, end) runs."""
    if len(arr) == 0:
        return []
    splits = np.where(np.diff(arr) > gap)[0] + 1
    return [(int(g[0]), int(g[-1])) for g in np.split(arr, splits)]


def neutral_dark_mask(img: np.ndarray, dark: int = 170, max_sat: int = 80) -> np.ndarray:
    """Fast uint8 mask of dark, near-neutral pixels (axes, spines, ticks, black curves)."""
    r, g, b = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    return (mx < dark) & ((mx - mn) < max_sat)


def find_box(
    img: np.ndarray,
    dark: int = 170,
    max_sat: int = 80,
    min_frac: float = 0.45,
    row_range: tuple[int, int] | None = None,
    col_range: tuple[int, int] | None = None,
) -> dict:
    """Locate the four plot-frame spines in a crop (optionally inside a sub-region)."""
    H, W = img.shape[:2]
    r0, r1 = row_range if row_range is not None else (0, H)
    c0, c1 = col_range if col_range is not None else (0, W)
    sub = img[r0:r1, c0:c1]
    dn = neutral_dark_mask(sub, dark=dark, max_sat=max_sat)
    rg = groups(np.where(dn.mean(axis=1) > min_frac)[0], gap=3)
    cg = groups(np.where(dn.mean(axis=0) > min_frac)[0], gap=3)
    if len(rg) < 2 or len(cg) < 2:
        raise RuntimeError(f"Could not find plot frame spines: rows={rg}, cols={cg}")
    top = (rg[0][0] + r0, rg[0][1] + r0)
    bot = (rg[-1][0] + r0, rg[-1][1] + r0)
    left = (cg[0][0] + c0, cg[0][1] + c0)
    right = (cg[-1][0] + c0, cg[-1][1] + c0)
    return {
        "top": top,
        "bot": bot,
        "left": left,
        "right": right,
        "r_top": (top[0] + top[1]) / 2.0 + 0.5,
        "r_bot": (bot[0] + bot[1]) / 2.0 + 0.5,
        "c_left": (left[0] + left[1]) / 2.0 + 0.5,
        "c_right": (right[0] + right[1]) / 2.0 + 0.5,
        "box": (float(left[1] + 1), float(top[1] + 1), float(right[0]), float(bot[0])),
    }


def _detect_ticks_1d(
    slice_2d: np.ndarray,
    min_len: int = 6,
    max_len: int = 250,
    max_width: int = 35,
) -> list[tuple[float, int, int]]:
    """Detect perpendicular tick runs starting within 3 px of a spine.
    slice_2d has shape (depth_from_spine, along_spine).
    Returns list of (center_offset, max_len, width).
    """
    depth, along = slice_2d.shape
    if depth < 3 or along == 0:
        return []
    lens = np.zeros(along, dtype=int)
    active = slice_2d[0] | slice_2d[1] | slice_2d[2]
    lens[active] = 1
    for k in range(1, depth):
        active = active & slice_2d[k]
        lens[active] += 1
    pos = np.where((lens >= min_len) & (lens <= max_len))[0]
    out = []
    for a, b in groups(pos, gap=2):
        w = b - a + 1
        if w <= max_width:
            out.append(((a + b) / 2.0 + 0.5, int(lens[a : b + 1].max()), w))
    return out


def detect_x_ticks(
    img: np.ndarray,
    frame: dict,
    side: str = "in",
    spine: str = "bottom",
    dark: int = 170,
    max_sat: int = 80,
    min_len: int = 6,
    max_len: int = 250,
    max_width: int = 35,
    depth: int = 120,
) -> list[tuple[float, int, int]]:
    """Detect x-axis ticks along the bottom or top spine, pointing in or out."""
    H, W = img.shape[:2]
    dn = neutral_dark_mask(img, dark=dark, max_sat=max_sat)
    c0, c1 = frame["left"][0], frame["right"][1] + 1
    if spine == "bottom":
        rb0, rb1 = frame["bot"]
        if side == "out":
            sl = dn[rb1 + 1 : min(H, rb1 + 1 + depth), c0:c1]
        else:
            sl = dn[max(0, rb0 - depth) : rb0, c0:c1][::-1]
    else:
        rt0, rt1 = frame["top"]
        if side == "out":
            sl = dn[max(0, rt0 - depth) : rt0, c0:c1][::-1]
        else:
            sl = dn[rt1 + 1 : min(H, rt1 + 1 + depth), c0:c1]
    return [(p + c0, l, w) for p, l, w in _detect_ticks_1d(sl, min_len, max_len, max_width)]


def detect_y_ticks(
    img: np.ndarray,
    frame: dict,
    side: str = "in",
    spine: str = "left",
    dark: int = 170,
    max_sat: int = 80,
    min_len: int = 6,
    max_len: int = 250,
    max_width: int = 35,
    depth: int = 120,
) -> list[tuple[float, int, int]]:
    """Detect y-axis ticks along the left or right spine, pointing in or out."""
    H, W = img.shape[:2]
    dn = neutral_dark_mask(img, dark=dark, max_sat=max_sat)
    r0, r1 = frame["top"][0], frame["bot"][1] + 1
    if spine == "left":
        cl0, cl1 = frame["left"]
        if side == "out":
            sl = dn[r0:r1, max(0, cl0 - depth) : cl0].T[::-1]
        else:
            sl = dn[r0:r1, cl1 + 1 : min(W, cl1 + 1 + depth)].T
    else:
        cr0, cr1 = frame["right"]
        if side == "out":
            sl = dn[r0:r1, cr1 + 1 : min(W, cr1 + 1 + depth)].T
        else:
            sl = dn[r0:r1, max(0, cr0 - depth) : cr0].T[::-1]
    return [(p + r0, l, w) for p, l, w in _detect_ticks_1d(sl, min_len, max_len, max_width)]


def select_regular_ticks(
    cands: list[tuple[float, int, int]],
    expected_n: int,
    lo_px: float | None = None,
    hi_px: float | None = None,
    tol_px: float = 8.0,
) -> list[float]:
    """Filter tick candidates by optional pixel range and pick the expected_n evenly spaced ticks."""
    pts = [p for p, _, _ in cands if (lo_px is None or p >= lo_px) and (hi_px is None or p <= hi_px)]
    pts = sorted(pts)
    if len(pts) == expected_n:
        return pts
    # Try every pair of endpoints (i, j) to find the best evenly spaced set of size expected_n
    best_set: list[float] | None = None
    best_err = 1e9
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            span = pts[j] - pts[i]
            if span < 50.0:
                continue
            step = span / (expected_n - 1)
            picked = []
            max_e = 0.0
            ok = True
            for k in range(expected_n):
                target = pts[i] + k * step
                near = min(pts, key=lambda p: abs(p - target))
                e = abs(near - target)
                if e > tol_px:
                    ok = False
                    break
                picked.append(near)
                max_e = max(max_e, e)
            if ok and len(set(picked)) == expected_n and max_e < best_err:
                best_err = max_e
                best_set = picked
    if best_set is None:
        raise RuntimeError(f"Could not select {expected_n} regular ticks from {pts}")
    return best_set


def calibrate_axis(positions: list[float], values: list[float], max_residual_px: float = 4.0) -> D.Axis:
    """Fit a linear Axis from detected tick pixel positions and printed labels."""
    if len(positions) != len(values):
        raise RuntimeError(f"Tick count mismatch: {len(positions)} positions vs {len(values)} labels")
    ax = D.Axis("linear", list(zip(positions, values)))
    if ax.residual_px > max_residual_px:
        raise RuntimeError(f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px")
    return ax


def col_runs(mask: np.ndarray, col: int, r_lo: int, r_hi: int) -> list[tuple[float, float]]:
    """Return list of (center_row, half_height) for contiguous True runs in column `col`."""
    if col < 0 or col >= mask.shape[1]:
        return []
    rr = np.where(mask[r_lo:r_hi, col])[0]
    if len(rr) == 0:
        return []
    out = []
    for a, b in groups(rr, gap=1):
        out.append((r_lo + (a + b) / 2.0 + 0.5, (b - a + 1) / 2.0))
    return out


def read_curve_at(
    mask: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    at: list[float],
    box: tuple[float, float, float, float],
    native_px: float = 1.0,
    search_cols: int = 4,
    band_fn: Callable[[float, float], tuple[float, float]] | None = None,
    pick: str = "closest_mid",
    note_prefix: str = "",
) -> list[dict]:
    """Read a solid curve from `mask` at requested x values `at`.
    `band_fn(x_val, col_px)` may return `(r_min, r_max)` to restrict the vertical search band.
    `pick` can be 'single', 'top', 'bottom', 'thickest', or 'closest_mid'.
    """
    c0, r0, c1, r1 = box
    r_lo_box = int(math.ceil(r0)) + 2
    r_hi_box = int(math.floor(r1)) - 2
    cal_x = max(xa.residual_px, 0.5 * native_px)
    cal_y = max(ya.residual_px, 0.5 * native_px)
    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        col = int(round(pc - 0.5))
        if band_fn is not None:
            b0, b1 = band_fn(xv, pc)
            r_lo = max(r_lo_box, int(math.floor(b0)))
            r_hi = min(r_hi_box, int(math.ceil(b1)))
        else:
            r_lo, r_hi = r_lo_box, r_hi_box

        runs, used = [], col
        for dc in sorted(range(-search_cols, search_cols + 1), key=abs):
            c_try = col + dc
            if c_try <= int(c0) + 1 or c_try >= int(c1) - 1:
                continue
            cand = col_runs(mask, c_try, r_lo, r_hi)
            if cand:
                runs, used = cand, c_try
                break

        if not runs:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"{note_prefix}no curve pixels at x={xv:g}"})
            continue

        if pick == "top":
            center, half = min(runs, key=lambda t: t[0])
        elif pick == "bottom":
            center, half = max(runs, key=lambda t: t[0])
        elif pick == "thickest":
            center, half = max(runs, key=lambda t: t[1])
        else:
            mid = 0.5 * (r_lo + r_hi)
            center, half = min(runs, key=lambda t: abs(t[0] - mid))

        # Local slope across +-step columns
        step = max(2, int(round(2 * native_px)))
        lc = col_runs(mask, used - step, max(r_lo, int(center - 60)), min(r_hi, int(center + 60)))
        rc = col_runs(mask, used + step, max(r_lo, int(center - 60)), min(r_hi, int(center + 60)))
        slope = 0.0
        if lc and rc:
            l_row = min(lc, key=lambda t: abs(t[0] - center))[0]
            r_row = min(rc, key=lambda t: abs(t[0] - center))[0]
            slope = (r_row - l_row) / (2.0 * step)
        offset = abs(pc - (used + 0.5))
        h_px = max(half, 0.5 * native_px) + 0.5 * native_px + abs(slope) * (offset + cal_x) + cal_y
        yv = ya.to_value(center)
        a, b = ya.to_value(center - h_px), ya.to_value(center + h_px)
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [float(pc), float(center)],
            "px_lo": float(center - h_px),
            "px_hi": float(center + h_px),
            "note": f"{note_prefix}cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px half={half:.1f}px",
        })
    return out


def read_dashed_at(
    mask: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    at: list[float],
    box: tuple[float, float, float, float],
    window_px: int = 160,
    sigma_px: float = 50.0,
    deg: int = 2,
    native_px: float = 1.0,
    band_fn: Callable[[float, float], tuple[float, float]] | None = None,
    min_run_h: float = 1.5,
    max_run_h: float = 40.0,
    note_prefix: str = "",
) -> list[dict]:
    """Read a dashed curve using local Gaussian-weighted polynomial regression across dashes.
    Avoids dash-gap misses and slanted dash-tip bias.
    """
    c0, r0, c1, r1 = box
    c_lo_box = int(math.ceil(c0)) + 2
    c_hi_box = int(math.floor(c1)) - 2
    r_lo_box = int(math.ceil(r0)) + 2
    r_hi_box = int(math.floor(r1)) - 2
    cal_x = max(xa.residual_px, 0.5 * native_px)
    cal_y = max(ya.residual_px, 0.5 * native_px)
    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        w_lo = max(c_lo_box, int(math.floor(pc - window_px)))
        w_hi = min(c_hi_box, int(math.ceil(pc + window_px)))
        xs_c, ys_r, halfs = [], [], []
        for c in range(w_lo, w_hi + 1):
            if band_fn is not None:
                b0, b1 = band_fn(xa.to_value(c + 0.5), c + 0.5)
                r_lo = max(r_lo_box, int(math.floor(b0)))
                r_hi = min(r_hi_box, int(math.ceil(b1)))
            else:
                r_lo, r_hi = r_lo_box, r_hi_box
            runs = [t for t in col_runs(mask, c, r_lo, r_hi) if min_run_h <= t[1] <= max_run_h]
            if len(runs) == 1:
                xs_c.append(c + 0.5)
                ys_r.append(runs[0][0])
                halfs.append(runs[0][1])
            elif len(runs) > 1:
                mid = 0.5 * (r_lo + r_hi)
                best = min(runs, key=lambda t: abs(t[0] - mid))
                xs_c.append(c + 0.5)
                ys_r.append(best[0])
                halfs.append(best[1])
        if len(xs_c) < deg + 2:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"{note_prefix}too few dash cols ({len(xs_c)}) near x={xv:g}"})
            continue
        dx = np.array(xs_c) - pc
        yr = np.array(ys_r)
        w = np.exp(-0.5 * (dx / sigma_px) ** 2)
        # Weighted polynomial fit (with iterative outlier trim for crossing/artifact immunity)
        use_deg = deg if (dx.max() - dx.min()) > 20.0 and len(dx) >= deg + 4 else min(1, deg)
        for _ in range(2):
            coeff = np.polyfit(dx, yr, use_deg, w=np.sqrt(w))
            resid = yr - np.polyval(coeff, dx)
            mad = float(np.median(np.abs(resid))) + 1.0
            keep = np.abs(resid) <= max(4.0 * mad, 8.0 * native_px)
            if keep.sum() >= use_deg + 2:
                dx, yr, w = dx[keep], yr[keep], w[keep]
        coeff = np.polyfit(dx, yr, use_deg, w=np.sqrt(w))
        center = float(coeff[-1])
        slope = float(coeff[-2]) if use_deg >= 1 else 0.0
        resid = yr - np.polyval(coeff, dx)
        wrms = float(math.sqrt(np.sum(w * resid * resid) / max(np.sum(w), 1e-12)))
        med_half = float(np.median(halfs))
        h_px = max(med_half, 0.5 * native_px) + max(wrms, 0.5 * native_px) + abs(slope) * cal_x + cal_y
        yv = ya.to_value(center)
        a, b = ya.to_value(center - h_px), ya.to_value(center + h_px)
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [float(pc), float(center)],
            "px_lo": float(center - h_px),
            "px_hi": float(center + h_px),
            "note": f"{note_prefix}dashed fit n={len(dx)} wrms={wrms:.1f}px cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


def find_errorbar_caps_at_x(
    color_mask: np.ndarray,
    xc: float,
    r_lo: int,
    r_hi: int,
    flank_inner: int = 30,
    flank_outer: int = 75,
    min_flank_px: int = 12,
    max_cap_h: int = 14,
) -> list[tuple[float, float, int]]:
    """Detect horizontal error bar caps at column `xc` by checking BOTH left and right flanks
    outside the circular marker radius. Returns list of `(center_row, half_height, width_sum)`.
    """
    H, W = color_mask.shape
    c = int(round(xc))
    l0, l1 = max(0, c - flank_outer), max(0, c - flank_inner)
    r0_c, r1_c = min(W, c + flank_inner), min(W, c + flank_outer)
    left_cnt = color_mask[r_lo:r_hi, l0:l1].sum(axis=1)
    right_cnt = color_mask[r_lo:r_hi, r0_c:r1_c].sum(axis=1)
    # Require strong horizontal cap presence outside marker radius (tolerant of a dashed line crossing one tip)
    is_cap = (
        ((left_cnt >= min_flank_px) & (right_cnt >= min_flank_px))
        | (((left_cnt + right_cnt) >= 2 * min_flank_px + 6) & ((left_cnt >= 4) | (right_cnt >= 4)))
    )
    rows = np.where(is_cap)[0]
    out = []
    for a, b in groups(rows, gap=2):
        h = b - a + 1
        if h <= max_cap_h:
            wsum = int((left_cnt[a : b + 1] + right_cnt[a : b + 1]).max())
            out.append((r_lo + (a + b) / 2.0 + 0.5, h / 2.0, wsum))
    return out
