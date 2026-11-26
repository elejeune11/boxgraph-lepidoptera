"""Shared image-based axis calibration and data extraction helpers for hansen2010 figures.

All functions operate purely on the passed RGB `image` numpy array (no PDF rendering or
external I/O). Axis calibrations are derived from detected frame lines and tick marks
directly in the crop.
"""
from __future__ import annotations

import math
from typing import Callable, Iterable

import numpy as np
from scipy import ndimage

from boxgraph import digitize as D


def runs_1d(arr: np.ndarray) -> list[tuple[int, int]]:
    """Inclusive (start, end) runs of True in a 1D boolean array."""
    idx = np.nonzero(arr)[0]
    if len(idx) == 0:
        return []
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    return [(int(s), int(e)) for s, e in zip(starts, ends)]


def gray(img: np.ndarray) -> np.ndarray:
    return img.astype(float).mean(axis=2)


def rgb_channels(img: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    f = img.astype(float)
    return f[:, :, 0], f[:, :, 1], f[:, :, 2]


def dark_mask(img: np.ndarray, thresh: float = 150.0) -> np.ndarray:
    return gray(img) < thresh


def black_ink_mask(img: np.ndarray, gray_max: float = 120.0, sat_max: float = 40.0) -> np.ndarray:
    r, g, b = rgb_channels(img)
    gr = (r + g + b) / 3.0
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    return (gr < gray_max) & (sat < sat_max)


def red_mask(img: np.ndarray, r_min: float = 120.0, rg_diff: float = 45.0, rb_diff: float = 45.0) -> np.ndarray:
    r, g, b = rgb_channels(img)
    return (r > r_min) & (r - g > rg_diff) & (r - b > rb_diff)


def blue_mask(img: np.ndarray, b_min: float = 110.0, br_diff: float = 35.0, bg_diff: float = 15.0) -> np.ndarray:
    r, g, b = rgb_channels(img)
    return (b > b_min) & (b - r > br_diff) & (b - g > bg_diff)


def green_mask(img: np.ndarray, g_min: float = 70.0, gr_diff: float = 25.0, gb_diff: float = 20.0) -> np.ndarray:
    r, g, b = rgb_channels(img)
    return (g > g_min) & (g - r > gr_diff) & (g - b > gb_diff)


def find_plot_frame(
    img: np.ndarray,
    thresh: float = 150.0,
    min_h_frac: float = 0.35,
    min_v_frac: float = 0.35,
    row_range: tuple[int, int] | None = None,
    col_range: tuple[int, int] | None = None,
) -> dict:
    """Locate plot frame spines from long dark horizontal and vertical runs.

    Returns dict with:
      cL, rT, cR, rB: inner plot edges (float)
      x_row, y_col: center of bottom x-axis and left y-axis spines
      top_row, right_col: center of top and right spines (if closed frame, else None)
      box: (cL, rT, cR, rB)
      x_span: (xc0, xc1), y_span: (yr0, yr1)
    """
    g = gray(img)
    H, W = g.shape
    r0_lim, r1_lim = row_range if row_range is not None else (0, H)
    c0_lim, c1_lim = col_range if col_range is not None else (0, W)
    d = np.zeros((H, W), dtype=bool)
    d[r0_lim:r1_lim, c0_lim:c1_lim] = g[r0_lim:r1_lim, c0_lim:c1_lim] < thresh

    w_span = c1_lim - c0_lim
    h_span = r1_lim - r0_lim

    min_h_px = min_h_frac * w_span
    min_v_px = min_v_frac * h_span
    row_sums = d[r0_lim:r1_lim, c0_lim:c1_lim].sum(axis=1)
    col_sums = d[r0_lim:r1_lim, c0_lim:c1_lim].sum(axis=0)

    long_rows = []
    for r_idx in np.nonzero(row_sums >= 0.85 * min_h_px)[0]:
        r = int(r0_lim + r_idx)
        runs = runs_1d(ndimage.binary_closing(d[r, c0_lim:c1_lim], structure=np.ones(5, bool)))
        if runs:
            s, e = max(runs, key=lambda x: x[1] - x[0])
            if e - s + 1 >= min_h_px:
                long_rows.append((r, e - s + 1, c0_lim + s, c0_lim + e))

    long_cols = []
    for c_idx in np.nonzero(col_sums >= 0.85 * min_v_px)[0]:
        c = int(c0_lim + c_idx)
        runs = runs_1d(ndimage.binary_closing(d[r0_lim:r1_lim, c], structure=np.ones(5, bool)))
        if runs:
            s, e = max(runs, key=lambda x: x[1] - x[0])
            if e - s + 1 >= min_v_px:
                long_cols.append((c, e - s + 1, r0_lim + s, r0_lim + e))

    if not long_rows or not long_cols:
        raise ValueError(f"Could not find plot frame spines (rows={len(long_rows)}, cols={len(long_cols)})")

    # Group contiguous spine bands
    def group_bands(items):
        items = sorted(items, key=lambda x: x[0])
        bands = []
        cur = [items[0]]
        for it in items[1:]:
            if it[0] <= cur[-1][0] + 6:
                cur.append(it)
            else:
                bands.append(cur)
                cur = [it]
        bands.append(cur)
        return bands

    rbands = group_bands(long_rows)
    cbands = group_bands(long_cols)

    # Bottom x-axis is lowest horizontal band; top spine (if present) is highest band
    bot_band = rbands[-1]
    top_band = rbands[0] if len(rbands) >= 2 else None
    left_band = cbands[0]
    right_band = cbands[-1] if len(cbands) >= 2 else None

    x_row = 0.5 * (bot_band[0][0] + bot_band[-1][0]) + 0.5
    y_col = 0.5 * (left_band[0][0] + left_band[-1][0]) + 0.5
    top_row = (0.5 * (top_band[0][0] + top_band[-1][0]) + 0.5) if top_band is not None else None
    right_col = (0.5 * (right_band[0][0] + right_band[-1][0]) + 0.5) if right_band is not None else None

    _, _, xc0, xc1 = max(bot_band, key=lambda x: x[1])
    _, _, yr0, yr1 = max(left_band, key=lambda x: x[1])

    cL = float(left_band[-1][0] + 1)
    rB = float(bot_band[0][0])
    rT = float(top_band[-1][0] + 1) if top_band is not None else float(yr0)
    cR = float(right_band[0][0]) if right_band is not None else float(xc1)

    return {
        "cL": cL,
        "rT": rT,
        "cR": cR,
        "rB": rB,
        "box": (cL, rT, cR, rB),
        "x_row": x_row,
        "y_col": y_col,
        "top_row": top_row,
        "right_col": right_col,
        "x_span": (xc0, xc1),
        "y_span": (yr0, yr1),
        "bot_band": (bot_band[0][0], bot_band[-1][0]),
        "left_band": (left_band[0][0], left_band[-1][0]),
        "top_band": (top_band[0][0], top_band[-1][0]) if top_band else None,
        "right_band": (right_band[0][0], right_band[-1][0]) if right_band else None,
    }


def detect_spine_ticks(
    img: np.ndarray,
    spine_pos: float,
    start: int,
    end: int,
    axis: str = "x",
    side: str = "out",
    thresh: float = 175.0,
    probe: tuple[int, int] = (3, 10),
    min_frac: float = 0.45,
    max_thick: int = 12,
    spine_band: tuple[int, int] | None = None,
) -> list[float]:
    """Detect tick mark centers along a horizontal ('x') or vertical ('y') axis spine."""
    g = gray(img)
    H, W = g.shape
    lo_off, hi_off = probe
    if axis == "x":
        s_lo, s_hi = spine_band if spine_band is not None else (int(round(spine_pos)), int(round(spine_pos)))
        if side == "out":
            rows = [r for r in range(s_hi + lo_off, s_hi + hi_off + 1) if 0 <= r < H]
        else:
            rows = [r for r in range(s_lo - hi_off, s_lo - lo_off + 1) if 0 <= r < H]
        start_c = max(0, int(start))
        end_c = min(W, int(end))
        if not rows or end_c <= start_c:
            return []
        prof = (g[rows, start_c:end_c] < thresh).mean(axis=0)
        offset0 = start_c
    else:
        s_lo, s_hi = spine_band if spine_band is not None else (int(round(spine_pos)), int(round(spine_pos)))
        if side == "out":
            cols = [c for c in range(s_lo - hi_off, s_lo - lo_off + 1) if 0 <= c < W]
        else:
            cols = [c for c in range(s_hi + lo_off, s_hi + hi_off + 1) if 0 <= c < W]
        start_r = max(0, int(start))
        end_r = min(H, int(end))
        if not cols or end_r <= start_r:
            return []
        prof = (g[start_r:end_r, cols] < thresh).mean(axis=1)
        offset0 = start_r

    ticks = []
    for s, e in runs_1d(prof >= min_frac):
        if e - s + 1 <= max_thick:
            # Intensity-weighted center for sub-pixel accuracy
            w = prof[s : e + 1]
            coords = np.arange(s, e + 1, dtype=float)
            center = offset0 + float((coords * w).sum() / w.sum()) + 0.5
            ticks.append(center)
    return ticks


def match_linear_axis(
    cands: list[float],
    values: list[float],
    increasing_px: bool = True,
    tol_px: float = 3.0,
    min_matches: int = 3,
    min_span_frac: float = 0.5,
    span_px: float | None = None,
    max_resid: float = 2.0,
    extra_anchors: list[float] | None = None,
) -> D.Axis:
    """Match detected tick candidates to known printed tick values on a linear axis.

    `values`: ordered printed tick values (increasing for x, increasing for y;
    for y-axis where pixel rows decrease as value increases, pass `increasing_px=False`).
    `extra_anchors`: optional spine corner positions (e.g., cL/cR or rB/rT) appended to
    candidates if a boundary tick coincides with a frame corner.
    """
    vals = sorted(float(v) for v in values)
    all_cands = sorted(set(float(c) for c in cands) | set(float(a) for a in (extra_anchors or []) if a is not None))
    if len(all_cands) < 2:
        raise D.CalibrationError(f"Too few tick candidates ({len(all_cands)}) for values {vals}")

    # Pairwise hypothesis search over (value_i, value_j) -> (cand_a, cand_b)
    v_span_total = vals[-1] - vals[0]
    best = None  # (n_matches, matched_v_span, -resid, marks)

    n_v = len(vals)
    n_c = len(all_cands)
    for i in range(n_v):
        for j in range(i + 1, n_v):
            dv = vals[j] - vals[i]
            if dv < 0.3 * v_span_total:
                continue
            for a in range(n_c):
                for b in range(n_c):
                    if a == b:
                        continue
                    dp = all_cands[b] - all_cands[a]
                    if increasing_px and dp <= 0:
                        continue
                    if not increasing_px and dp >= 0:
                        continue
                    slope = dp / dv  # px per unit value
                    total_px = abs(slope * v_span_total)
                    if span_px is not None:
                        if total_px < min_span_frac * span_px or total_px > 1.15 * span_px:
                            continue
                    # Match all values
                    marks = []
                    for v in vals:
                        pred_p = all_cands[a] + slope * (v - vals[i])
                        nearest = min(all_cands, key=lambda c: abs(c - pred_p))
                        if abs(nearest - pred_p) <= tol_px:
                            marks.append((nearest, v))
                    if len(marks) < min(min_matches, n_v):
                        continue
                    # Check uniqueness of pixels
                    if len({round(m[0], 3) for m in marks}) < len(marks):
                        continue
                    ax = D.Axis("linear", marks)
                    if ax.residual_px > max_resid:
                        continue
                    v_cov = marks[-1][1] - marks[0][1]
                    score = (len(marks), round(v_cov, 6), -ax.residual_px)
                    if best is None or score > best[0]:
                        best = (score, ax, marks)

    if best is None:
        raise D.CalibrationError(
            f"Could not match linear axis for values {vals} from {len(all_cands)} candidates: "
            f"{[round(c, 1) for c in all_cands]}"
        )
    return best[1]


def val_interval(y_axis: D.Axis, row: float, half_px: float) -> tuple[float, float, float]:
    """Return (y, lo, hi) at pixel `row` with symmetric pixel uncertainty `half_px`."""
    y = y_axis.to_value(row)
    a = y_axis.to_value(row - half_px)
    b = y_axis.to_value(row + half_px)
    lo, hi = min(a, b), max(a, b)
    lo = min(lo, y)
    hi = max(hi, y)
    return float(y), float(lo), float(hi)


def sample_curve_column(
    mask: np.ndarray,
    col_f: float,
    r_min: int,
    r_max: int,
    half_w: int = 2,
    band: tuple[int, int] | None = None,
) -> tuple[float, float, float, int] | None:
    """Sample a binary curve mask around float column `col_f`.

    Returns (center_row, top_row, bot_row, count) or None if no ink found.
    If `band` = (r_lo, r_hi) is given, restricts vertical search to [r_lo, r_hi].
    """
    H, W = mask.shape
    c0 = max(0, int(round(col_f)) - half_w)
    c1 = min(W, int(round(col_f)) + half_w + 1)
    rr0 = max(0, r_min if band is None else max(r_min, band[0]))
    rr1 = min(H, r_max if band is None else min(r_max, band[1]))
    if rr1 <= rr0 or c1 <= c0:
        return None
    sub = mask[rr0:rr1, c0:c1]
    ys = np.nonzero(sub)[0]
    if len(ys) == 0:
        return None
    ys_abs = ys + rr0
    return (
        float(np.median(ys_abs)) + 0.5,
        float(ys_abs.min()),
        float(ys_abs.max() + 1),
        int(len(ys_abs)),
    )


def trace_curve_at(
    mask: np.ndarray,
    x_axis: D.Axis,
    y_axis: D.Axis,
    box: tuple[float, float, float, float],
    at: Iterable[float],
    half_w: int = 2,
    inset_px: int = 3,
    extra_err_px: float = 1.0,
    search_w_fallback: int = 8,
    band_fn: Callable[[float], tuple[int, int]] | None = None,
    select: str = "median",
) -> list[dict]:
    """Read a single-color curve mask at every x in `at`."""
    cL, rT, cR, rB = box
    r_min = int(math.ceil(rT + inset_px))
    r_max = int(math.floor(rB - inset_px))
    cal_err = max(y_axis.residual_px, 0.5)
    out = []
    for x in at:
        x_val = float(x)
        c_f = x_axis.to_pixel(x_val)
        band = band_fn(x_val) if band_fn is not None else None
        # Try narrow window first, then fallback slightly wider for dashed/dotted curves
        hit = None
        used_w = half_w
        for hw in (half_w, search_w_fallback):
            c0 = max(int(math.ceil(cL + 1)), int(round(c_f)) - hw)
            c1 = min(int(math.floor(cR - 1)), int(round(c_f)) + hw + 1)
            rr0 = max(0, r_min if band is None else max(r_min, band[0]))
            rr1 = min(mask.shape[0], r_max if band is None else min(r_max, band[1]))
            if rr1 <= rr0 or c1 <= c0:
                continue
            sub = mask[rr0:rr1, c0:c1]
            ys = np.nonzero(sub)[0]
            if len(ys) > 0:
                ys_abs = ys + rr0
                if select == "median":
                    r_c = float(np.median(ys_abs)) + 0.5
                    r_lo_px = float(np.percentile(ys_abs, 10))
                    r_hi_px = float(np.percentile(ys_abs, 90)) + 1.0
                elif select == "mean":
                    r_c = float(np.mean(ys_abs)) + 0.5
                    r_lo_px = float(ys_abs.min())
                    r_hi_px = float(ys_abs.max()) + 1.0
                elif select == "upper":
                    # Topmost contiguous vertical run (smallest row)
                    r_runs = runs_1d(sub.any(axis=1))
                    s, e = r_runs[0]
                    r_c = rr0 + 0.5 * (s + e) + 0.5
                    r_lo_px = float(rr0 + s)
                    r_hi_px = float(rr0 + e + 1)
                elif select == "lower":
                    r_runs = runs_1d(sub.any(axis=1))
                    s, e = r_runs[-1]
                    r_c = rr0 + 0.5 * (s + e) + 0.5
                    r_lo_px = float(rr0 + s)
                    r_hi_px = float(rr0 + e + 1)
                else:
                    raise ValueError(f"unknown select={select!r}")
                hit = (r_c, r_lo_px, r_hi_px, int(len(ys_abs)))
                used_w = hw
                break
        if hit is None:
            out.append({"x": x_val, "y": None, "note": f"no curve pixels near col {c_f:.1f}"})
            continue
        r_c, r_top, r_bot, cnt = hit
        half_thick = 0.5 * max(r_bot - r_top, 1.0)
        # Include local slope contribution if fallback widened x
        extra = extra_err_px + (0.3 * used_w if used_w > half_w else 0.0)
        half_px = half_thick + cal_err + extra
        y, lo, hi = val_interval(y_axis, r_c, half_px)
        out.append({
            "x": x_val,
            "y": y,
            "lo": lo,
            "hi": hi,
            "px": [float(c_f), float(r_c)],
            "px_lo": float(r_c - half_px),
            "px_hi": float(r_c + half_px),
            "note": f"curve {cnt}px (thick {r_bot - r_top:.1f}px, cal {cal_err:.2f}px)",
        })
    return out


def measure_error_bar_vertical(
    img: np.ndarray,
    col_c: float,
    r_center: float,
    r_search_top: int,
    r_search_bot: int,
    y_axis: D.Axis,
    dark_thresh: float = 125.0,
    stem_half_w: int = 3,
    gap_tol: int = 3,
    mirror_lower_if_hidden: bool = True,
    bar_fill_gray: float | None = None,
) -> dict:
    """Trace a vertical black error bar upward and downward from (col_c, r_center).

    If `bar_fill_gray` is given (e.g. a gray bar with a darker black error bar drawn inside),
    uses a stricter threshold below the bar top so the gray fill is not mistaken for stem.
    """
    g = gray(img)
    H, W = g.shape
    c0 = max(0, int(round(col_c)) - stem_half_w)
    c1 = min(W, int(round(col_c)) + stem_half_w + 1)

    # Trace UP from r_center - 1
    up_thresh = dark_thresh
    r_up = int(math.floor(r_center))
    best_up = r_up
    gap = 0
    for r in range(r_up - 1, max(0, r_search_top) - 1, -1):
        # Check if dark ink exists near center column
        if (g[r, c0:c1] < up_thresh).any():
            best_up = r
            gap = 0
        else:
            gap += 1
            if gap > gap_tol:
                break

    # Refine center of horizontal top cap (if cap is 2-4 px thick)
    # Find run of rows near best_up that are wider than stem (cap) or just use best_up + 0.5
    r_cap_top = float(best_up) + 0.5
    # Check if there is a slightly thicker dark cap right at best_up..best_up+3
    cap_rows = [r for r in range(best_up, min(r_up, best_up + 4)) if (g[r, max(0, c0 - 3) : min(W, c1 + 3)] < up_thresh).sum() >= 5]
    if cap_rows:
        r_cap_top = 0.5 * (min(cap_rows) + max(cap_rows)) + 0.5

    # Trace DOWN from r_center + 1
    down_thresh = min(dark_thresh, bar_fill_gray - 35.0) if bar_fill_gray is not None else dark_thresh
    r_dn = int(math.ceil(r_center))
    best_dn = r_dn
    found_down = False
    if down_thresh > 25.0:
        gap = 0
        # Start slightly below top edge border of bar (3 px) if bar_fill_gray is set
        start_dn = r_dn + (4 if bar_fill_gray is not None else 1)
        for r in range(start_dn, min(H, r_search_bot)):
            if (g[r, c0:c1] < down_thresh).any():
                best_dn = r
                found_down = True
                gap = 0
            else:
                gap += 1
                if gap > gap_tol:
                    break

    cal_err = max(y_axis.residual_px, 0.5)
    err_px = 1.5 + cal_err
    y_val = y_axis.to_value(r_center)

    bar_hi = float(y_axis.to_value(r_cap_top))
    bar_hi_outer = float(y_axis.to_value(r_cap_top - err_px))

    if found_down and best_dn > r_dn + 4:
        cap_dn_rows = [r for r in range(max(r_dn, best_dn - 3), best_dn + 1) if (g[r, max(0, c0 - 3) : min(W, c1 + 3)] < down_thresh).sum() >= 5]
        r_cap_bot = (0.5 * (min(cap_dn_rows) + max(cap_dn_rows)) + 0.5) if cap_dn_rows else (float(best_dn) + 0.5)
        bar_lo = float(y_axis.to_value(r_cap_bot))
        bar_lo_outer = float(y_axis.to_value(r_cap_bot + err_px))
        note = f"err bar [{r_cap_top:.1f}, {r_cap_bot:.1f}]px"
    elif mirror_lower_if_hidden:
        bar_lo = float(2.0 * y_val - bar_hi)
        bar_lo_outer = float(2.0 * y_val - bar_hi_outer)
        r_cap_bot = float(2.0 * r_center - r_cap_top)
        note = f"upper err bar {r_cap_top:.1f}px (lower hidden inside dark bar: mirrored)"
    else:
        bar_lo = float(y_val)
        bar_lo_outer = float(y_axis.to_value(r_center + err_px))
        r_cap_bot = float(r_center)
        note = f"upper err bar {r_cap_top:.1f}px (no lower bar found)"

    return {
        "bar_lo": min(bar_lo, bar_hi),
        "bar_hi": max(bar_lo, bar_hi),
        "bar_lo_outer": min(bar_lo_outer, bar_hi_outer),
        "bar_hi_outer": max(bar_lo_outer, bar_hi_outer),
        "px_lo": float(r_cap_top),
        "px_hi": float(r_cap_bot),
        "note": note,
    }


# ==========================================================================
# Batch 1 helpers: Figs. 1C, 1D, 2D, 2E, 2H
# ==========================================================================


def calibrate_fig1c(img: np.ndarray):
    fr = find_plot_frame(img)
    xt = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "in", spine_band=fr["bot_band"])
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    xa = match_linear_axis(xt, [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600], increasing_px=True, extra_anchors=[fr["y_col"]])
    ya = match_linear_axis(yt, [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600], increasing_px=False, extra_anchors=[fr["x_row"]])
    return xa, ya, fr


def read_fig1c_series(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, fr = calibrate_fig1c(img)
    r, g, b = rgb_channels(img)
    if series == "actin_alone":
        mask = (r < 125) & (g < 125) & (b < 125) & (np.abs(r - b) < 35)
    elif series == "vasp":
        mask = (r > 125) & (r - b > 30) & (r - g > 20)
    elif series == "kck_vasp":
        mask = (b > 115) & (b - r > 25) & (b - g > 15)
    else:
        raise ValueError(f"unknown Fig. 1C series {series!r}")
    # Exclude text labels above the curves: "KCK-VASP" / "VASP" around x=300..900, y>950; "Actin alone" around x=700..1300, y=500..750
    cL, rT, cR, rB = fr["box"]
    if series == "kck_vasp":
        mask[: int(ya.to_pixel(1150)), int(xa.to_pixel(250)) : int(xa.to_pixel(1050))] = False
    elif series == "vasp":
        mask[: int(ya.to_pixel(1000)), int(xa.to_pixel(250)) : int(xa.to_pixel(850))] = False
    elif series == "actin_alone":
        mask[: int(ya.to_pixel(450)), int(xa.to_pixel(600)) : int(xa.to_pixel(1400))] = False
    return trace_curve_at(mask, xa, ya, (cL, rT, cR, rB), at, half_w=5, inset_px=2, extra_err_px=1.5)


def detect_fig1d_bars(img: np.ndarray):
    """Detect all 14 vertical bars in Fig. 1D and calibrate y-axis (0..35 sub/s)."""
    fr = find_plot_frame(img)
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(
        yt,
        [0, 5, 10, 15, 20, 25, 30, 35],
        increasing_px=False,
        extra_anchors=[fr["x_row"], fr["top_row"]],
    )
    g = gray(img)
    rB = int(round(fr["rB"]))
    rT = int(round(fr["rT"]))
    cL = int(round(fr["cL"])) + 4
    cR = int(round(fr["cR"])) - 4

    # 1. White & mid-gray bars have two vertical black outlines standing on rB
    v_dark = (g[rB - 90 : rB - 10, cL:cR] < 130).mean(axis=0) > 0.75
    dark_runs = [(cL + s, cL + e) for s, e in runs_1d(v_dark)]

    bars = []  # (c_left, c_right, kind: 'white'|'gray'|'dark')
    # Solid dark charcoal bars (+Cy3-VASP) appear as wide dark runs (50..80 px)
    thin_cols = []
    for s, e in dark_runs:
        w = e - s + 1
        if 45 <= w <= 85:
            bars.append((s, e, "dark"))
        elif w <= 12:
            thin_cols.append(0.5 * (s + e))

    # Pair thin vertical borders separated by 50..75 px
    i = 0
    while i + 1 < len(thin_cols):
        w = thin_cols[i + 1] - thin_cols[i]
        if 48 <= w <= 78:
            c0, c1 = int(round(thin_cols[i])), int(round(thin_cols[i + 1]))
            mg = float(np.median(g[rB - 80 : rB - 15, c0 + 10 : c1 - 10]))
            kind = "white" if mg > 210 else "gray"
            bars.append((c0, c1, kind))
            i += 2
        else:
            i += 1

    bars.sort(key=lambda x: x[0])
    if len(bars) != 14:
        raise ValueError(f"Fig. 1D: expected 14 bars, found {len(bars)}: {bars}")

    cal_err = max(ya.residual_px, 0.5)
    measured = []
    for c0, c1, kind in bars:
        c_mid = 0.5 * (c0 + c1) + 0.5
        # Find top edge of bar by scanning down from top in left/right quarters of bar (avoiding center error-bar stem!)
        q_cols = list(range(c0 + 8, int(c_mid) - 10)) + list(range(int(c_mid) + 10, c1 - 7))
        # Scan upward from rB - 20 until bar interior ends (white above bar)
        # For white bar, top edge is a horizontal black border line spanning q_cols!
        if kind == "white":
            # Look for horizontal dark border line above rB - 30 spanning q_cols
            prof = (g[rT + 20 : rB - 30, q_cols] < 130).mean(axis=1)
            r_hits = runs_1d(prof >= 0.7)
            if not r_hits:
                raise ValueError(f"Fig. 1D: could not locate top border of white bar [{c0},{c1}]")
            s, e = r_hits[-1]  # lowest horizontal border above baseline (below legend!)
            r_top = (rT + 20) + 0.5 * (s + e) + 0.5
            fill_gray = 250.0
        else:
            # Gray or dark bar: find highest row where fill+border is darker than 200 across >= 75% of q_cols
            # Start from bottom up, finding where it transitions to white (>230)
            r_cur = rB - 20
            while r_cur > rT + 20:
                if (g[r_cur, q_cols] < 190).mean() < 0.5:
                    break
                r_cur -= 1
            # r_cur+1 is top edge of bar
            r_top = float(r_cur + 1) + 0.5
            fill_gray = float(np.median(g[int(r_top) + 15 : rB - 15, q_cols]))

        y, lo, hi = val_interval(ya, r_top, 1.5 + cal_err)
        eb = measure_error_bar_vertical(
            img,
            c_mid,
            r_top,
            rT + 4,
            rB - 4,
            ya,
            dark_thresh=125.0,
            stem_half_w=4,
            gap_tol=3,
            mirror_lower_if_hidden=True,
            bar_fill_gray=fill_gray if kind != "white" else None,
        )
        measured.append({
            "kind": kind,
            "c_mid": c_mid,
            "r_top": r_top,
            "y": y,
            "lo": lo,
            "hi": hi,
            **eb,
        })
    return measured, ya, fr


def read_fig1d_series(img: np.ndarray, series: str) -> list[dict]:
    bars, _, _ = detect_fig1d_bars(img)
    # Bars 0..10 are 50 mM KCl:
    #   0: 0 nM (white)
    #   1,2: 1 nM (gray=+VASP, dark=+Cy3-VASP)
    #   3,4: 5 nM
    #   5,6: 10 nM
    #   7,8: 50 nM
    #   9,10: 100 nM
    # Bars 11..13 are 100 mM KCl:
    #   11: 0 nM (white), 12: 50 nM (gray), 13: 100 nM (gray)
    if series == "50mM_vasp":
        picks = [(0.0, 0), (1.0, 1), (5.0, 3), (10.0, 5), (50.0, 7), (100.0, 9)]
    elif series == "50mM_cy3_vasp":
        picks = [(1.0, 2), (5.0, 4), (10.0, 6), (50.0, 8), (100.0, 10)]
    elif series == "100mM_vasp":
        picks = [(0.0, 11), (50.0, 12), (100.0, 13)]
    else:
        raise ValueError(f"unknown Fig. 1D series {series!r}")
    out = []
    for conc, idx in picks:
        b = bars[idx]
        out.append({
            "x": conc,
            "y": b["y"],
            "lo": b["lo"],
            "hi": b["hi"],
            "bar_lo": b["bar_lo"],
            "bar_hi": b["bar_hi"],
            "bar_lo_outer": b["bar_lo_outer"],
            "bar_hi_outer": b["bar_hi_outer"],
            "px": [b["c_mid"], b["r_top"]],
            "px_lo": b["px_lo"],
            "px_hi": b["px_hi"],
            "note": f"{b['kind']} bar; {b['note']}",
        })
    return out


def calibrate_fig2d(img: np.ndarray):
    fr = find_plot_frame(img)
    xt = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "in", spine_band=fr["bot_band"])
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    xa = match_linear_axis(xt, [0, 10, 20, 30, 40, 50], increasing_px=True, extra_anchors=[fr["right_col"]])
    ya = match_linear_axis(yt, [0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False, extra_anchors=[fr["top_row"]])
    return xa, ya, fr


def calibrate_fig2e_inset(img: np.ndarray):
    fr = find_plot_frame(img)
    xt = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "in", spine_band=fr["bot_band"])
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    xa = match_linear_axis(xt, [0, 1, 2, 3, 4, 5], increasing_px=True)
    ya = match_linear_axis(yt, [0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False)
    return xa, ya, fr


def calibrate_fig2h(img: np.ndarray):
    fr = find_plot_frame(img)
    xt = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "in", spine_band=fr["bot_band"])
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    xa = match_linear_axis(xt, [0, 0.5, 1.0, 1.5, 2.0, 2.5], increasing_px=True, extra_anchors=[fr["right_col"]])
    ya = match_linear_axis(yt, [0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False, extra_anchors=[fr["top_row"]])
    return xa, ya, fr


# ==========================================================================
# Batch 2 helpers: Figs. 3C, 3D, 3F, 3G, 3I, 4E, 4G, 4H, 4I
# ==========================================================================


def calibrate_fig3c(img: np.ndarray):
    # Exclude horizontal bracket line below the x-axis (located > 72% down crop)
    H = img.shape[0]
    fr = find_plot_frame(img, row_range=(0, int(0.70 * H)))
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yt, [0.0, 0.05, 0.10, 0.15, 0.20, 0.25], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 1.0), (fr["right_col"], 10.0)])
    return xa, ya, fr


def read_fig3c_series(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, fr = calibrate_fig3c(img)
    r, g, b = rgb_channels(img)
    if series == "v":
        mask = (r < 90) & (g < 90) & (b < 90) & (np.abs(r - g) < 25) & (np.abs(r - b) < 25)
    elif series == "v_cp":
        mask = (g > 75) & (g - r > 20) & (g - b > 15)
    elif series == "v_a":
        mask = (r > 130) & (r - g > 38) & (r - b > 38)
    elif series == "v_a_cp":
        mask = (b > 105) & (b - r > 25) & (b - g > 15)
    else:
        raise ValueError(f"unknown Fig. 3C series {series!r}")
    # Exclude legend in upper right (x > 4.5, y > 0.16)
    mask[: int(ya.to_pixel(0.16)), int(xa.to_pixel(4.5)) :] = False
    cL, rT, cR, rB = fr["box"]
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for x in at:
        xv = float(x)
        c_f = xa.to_pixel(xv)
        # At pixel 1 (left spine) or 10 (right spine), sample narrowly right inside the spine
        if xv <= 1.01:
            c0, c1 = int(round(cL)) + (3 if series == "v" else 1), int(round(cL)) + 7
        elif xv >= 9.99:
            c0, c1 = int(round(cR)) - 7, int(round(cR)) - 1
        else:
            c0, c1 = max(int(cL) + 1, int(round(c_f)) - 4), min(int(cR) - 1, int(round(c_f)) + 5)
        sub = mask[int(rT) + 2 : int(rB) - 2, c0:c1]
        ys = np.nonzero(sub)[0]
        if len(ys) == 0:
            # Widen slightly if another curve sits right on top of this color at a crossing
            c0w, c1w = max(int(cL) + 1, int(round(c_f)) - 32), min(int(cR) - 1, int(round(c_f)) + 33)
            sub = mask[int(rT) + 2 : int(rB) - 2, c0w:c1w]
            ys = np.nonzero(sub)[0]
        if len(ys) == 0:
            out.append({"x": xv, "y": None, "note": f"no {series} pixels near pixel {xv}"})
            continue
        ys_abs = int(rT) + 2 + ys
        r_c = float(np.median(ys_abs)) + 0.5
        thick = float(np.percentile(ys_abs, 90) - np.percentile(ys_abs, 10)) + 1.0
        half_px = 0.5 * thick + cal_err + 1.5
        y, lo, hi = val_interval(ya, r_c, half_px)
        out.append({
            "x": xv, "y": y, "lo": lo, "hi": hi,
            "px": [float(c_f), float(r_c)],
            "px_lo": float(r_c - half_px), "px_hi": float(r_c + half_px),
            "note": f"Fig. 3C {series} ({len(ys)} px)",
        })
    return out


def read_fig3d_group(img: np.ndarray, group: str) -> list[dict]:
    """Read Barbed end (bars 0,1) or Lateral side (bars 2,3) dwell times with SEM error bars."""
    fr = find_plot_frame(img)
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yt, [0.0, 0.2, 0.4, 0.6, 0.8], increasing_px=False, extra_anchors=[fr["x_row"]])
    g = gray(img)
    rB, rT = int(round(fr["rB"])), int(round(fr["rT"]))
    cL, cR = int(round(fr["cL"])) + 4, int(round(fr["cR"])) - 4
    non_w = (g[rB - 35 : rB - 8, cL:cR] < 235).mean(axis=0) > 0.7
    bcols = [(cL + s, cL + e) for s, e in runs_1d(non_w) if e - s >= 40]
    if len(bcols) != 4:
        raise ValueError(f"Fig. 3D: expected 4 bars, found {len(bcols)}: {bcols}")
    idxs = [0, 1] if group == "barbed_end" else [2, 3]
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for k, bi in enumerate(idxs, start=1):
        c0, c1 = bcols[bi]
        c_mid = 0.5 * (c0 + c1) + 0.5
        q = list(range(c0 + 10, int(c_mid) - 12)) + list(range(int(c_mid) + 12, c1 - 10))
        # Scan upward from bottom until fill turns white (>230)
        r_cur = rB - 8
        while r_cur > rT + 10:
            if (g[r_cur, q] < 225).mean() < 0.5:
                break
            r_cur -= 1
        r_top = float(r_cur + 1) + 0.5
        fill_g = float(np.median(g[int(r_top) + 10 : rB - 5, q]))
        y, lo, hi = val_interval(ya, r_top, 1.5 + cal_err)
        eb = measure_error_bar_vertical(
            img, c_mid, r_top, rT + 4, rB - 2, ya,
            dark_thresh=115.0, stem_half_w=4, gap_tol=3,
            mirror_lower_if_hidden=True, bar_fill_gray=fill_g,
        )
        cond_lbl = "+Actin" if k == 1 else "+Actin, CP"
        out.append({
            "x": float(k), "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_mid, r_top], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"{group} {cond_lbl}; {eb['note']}",
        })
    return out


def calibrate_fig3f(img: np.ndarray):
    fr = find_plot_frame(img)
    yt_l = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    yt_r = detect_spine_ticks(img, fr["right_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "out", spine_band=fr["right_band"])
    ya_l = match_linear_axis(yt_l, [0, 5, 10, 15, 20, 25], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    ya_r = match_linear_axis(yt_r, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    # Calibrate x from the 9 printed numeric tick labels (4, 6, 8, 10, 12, 14, 16, 18, 20) directly below x_row
    g = gray(img)
    rB = int(round(fr["x_row"]))
    c0, c1 = max(0, int(fr["y_col"]) - 30), min(img.shape[1], int(fr["right_col"]) + 30)
    sub = g[rB + 12 : min(img.shape[0], rB + 58), c0:c1] < 160
    col_has = ndimage.binary_closing(sub.any(axis=0), structure=np.ones(16, bool))
    lbls = [c0 + 0.5 * (s + e) + 0.5 for s, e in runs_1d(col_has) if e - s >= 5]
    xa = match_linear_axis(lbls, [4, 6, 8, 10, 12, 14, 16, 18, 20], increasing_px=True, tol_px=8.0, max_resid=4.5)
    return xa, ya_l, ya_r, fr


def calibrate_fig3g(img: np.ndarray):
    fr = find_plot_frame(img)
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yt, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 8.0)])
    return xa, ya, fr


def calibrate_fig3i(img: np.ndarray):
    fr = find_plot_frame(img)
    yt_r = detect_spine_ticks(img, fr["right_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "out", spine_band=fr["right_band"])
    ya_r = match_linear_axis(yt_r, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 2.0), (fr["right_col"], 18.0)])
    return xa, ya_r, fr


def calibrate_fig4e_inset(img: np.ndarray):
    fr = find_plot_frame(img)
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yt, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False, extra_anchors=[fr["x_row"]])
    # On Fig. 4E inset, x = 0 sits at y_col and x = 20 tick sits slightly left of right_col (detected via inward tick/labels)
    xt = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "in", thresh=185, probe=(2, 8), spine_band=fr["bot_band"])
    # Also locate x-axis tick labels (0, 5, 10, 15, 20) below x_row
    g = gray(img)
    rB = int(round(fr["x_row"]))
    c0, c1 = max(0, int(fr["y_col"]) - 30), min(img.shape[1], int(fr["right_col"]) + 20)
    sub = g[rB + 12 : min(img.shape[0], rB + 58), c0:c1] < 160
    col_has = ndimage.binary_closing(sub.any(axis=0), structure=np.ones(20, bool))
    lbls = [c0 + 0.5 * (s + e) + 0.5 for s, e in runs_1d(col_has) if e - s >= 5]
    # Anchor 0 at y_col and 20 at detected x=20 tick (~1111) or label
    t20 = [c for c in xt if fr["y_col"] + 0.85 * (fr["right_col"] - fr["y_col"]) < c < fr["right_col"] - 15]
    cands = [fr["y_col"]] + t20 + lbls
    xa = match_linear_axis(cands, [0, 5, 10, 15, 20], increasing_px=True, tol_px=14.0, max_resid=6.0)
    return xa, ya, fr


def calibrate_fig4g(img: np.ndarray):
    fr = find_plot_frame(img)
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yt, [5, 10, 15, 20, 25, 30, 35], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 120.0)])
    return xa, ya, fr


def read_fig4h_series(img: np.ndarray, series: str) -> list[dict]:
    """Read overlapping bars in Fig. 4H: dark gray 'no VASP' and light gray '+ 50 nM VASP'."""
    fr = find_plot_frame(img)
    yt = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yt, [0, 10, 20, 30, 40, 50, 60], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    g = gray(img)
    rB, rT = int(round(fr["rB"])), int(round(fr["rT"]))
    cL, cR = int(round(fr["cL"])) + 4, int(round(fr["cR"])) - 4
    non_w = (g[rB - 55 : rB - 8, cL:cR] < 235).mean(axis=0) > 0.7
    bcols = [(cL + s, cL + e) for s, e in runs_1d(non_w) if e - s >= 45]
    concs = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
    if len(bcols) != len(concs):
        raise ValueError(f"Fig. 4H: expected 6 bar groups, found {len(bcols)}: {bcols}")
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for conc, (c0, c1) in zip(concs, bcols):
        c_mid = 0.5 * (c0 + c1) + 0.5
        q = list(range(c0 + 10, int(c_mid) - 12)) + list(range(int(c_mid) + 12, c1 - 10))
        prof = np.median(g[rT + 10 : rB - 2, q], axis=1)
        if series == "no_vasp":
            # Bottommost contiguous dark-gray run (< 165) standing on rB
            runs = runs_1d(prof < 165.0)
            r_top = (rT + 10) + float(runs[-1][0]) + 0.5
            fill_g = float(np.median(g[int(r_top) + 4 : rB - 2, q])) if int(r_top) + 6 < rB else 135.0
            # At 0.25 uM the inner bar is tiny (~1.1 sub/s = 16 px)
            eb_thresh = 100.0
        elif series == "50nM_vasp":
            # Bottommost contiguous non-white run (< 225) standing on rB (ignores legend above bars 0..3!)
            runs = runs_1d(prof < 225.0)
            r_top = (rT + 10) + float(runs[-1][0]) + 0.5
            fill_g = 205.0
            eb_thresh = 120.0
        else:
            raise ValueError(f"unknown Fig. 4H series {series!r}")
        y, lo, hi = val_interval(ya, r_top, 1.5 + cal_err)
        eb = measure_error_bar_vertical(
            img, c_mid, r_top, rT + 4, rB - 2, ya,
            dark_thresh=eb_thresh, stem_half_w=4, gap_tol=3,
            mirror_lower_if_hidden=True, bar_fill_gray=fill_g,
        )
        out.append({
            "x": conc, "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_mid, r_top], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"Fig. 4H {series}; {eb['note']}",
        })
    return out


def calibrate_fig4i(img: np.ndarray):
    fr = find_plot_frame(img)
    xt = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "in", spine_band=fr["bot_band"])
    yt_l = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    yt_r = detect_spine_ticks(img, fr["right_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "out", spine_band=fr["right_band"])
    xa = match_linear_axis(xt, [0, 0.5, 1.0, 1.5, 2.0, 2.5], increasing_px=True, extra_anchors=[fr["y_col"], fr["right_col"]])
    ya_l = match_linear_axis(yt_l, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    ya_r = match_linear_axis(yt_r, [0, 10, 20, 30, 40, 50, 60, 70], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    return xa, ya_l, ya_r, fr


# ==========================================================================
# Batch 3 helpers: Figs. 5A, 5B, 5C, 5D, 5E
# ==========================================================================


def calibrate_fig5a(img: np.ndarray):
    # Top absorbance panel ends around 80% down the crop (above the residuals panel)
    H = img.shape[0]
    fr = find_plot_frame(img, row_range=(0, int(0.81 * H)))
    xo = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "out", spine_band=fr["bot_band"])
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    xa = match_linear_axis(xo, [5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0], increasing_px=True, extra_anchors=[fr["y_col"], fr["right_col"]])
    ya = match_linear_axis(yi, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False, extra_anchors=[fr["top_row"]])
    return xa, ya, fr


def read_fig5a_curve(img: np.ndarray, at: list[float], color_name: str) -> list[dict]:
    xa, ya, fr = calibrate_fig5a(img)
    if color_name == "green":
        mask = green_mask(img, g_min=75.0, gr_diff=18.0, gb_diff=15.0)
    elif color_name == "blue":
        mask = blue_mask(img, b_min=105.0, br_diff=25.0, bg_diff=12.0)
    elif color_name == "red":
        mask = red_mask(img, r_min=120.0, rg_diff=35.0, rb_diff=35.0)
    else:
        raise ValueError(f"unknown Fig. 5A color {color_name!r}")
    # Exclude cartoons at the top (y > 0.85)
    mask[: int(ya.to_pixel(0.85)), :] = False
    cL, rT, cR, rB = fr["box"]
    rB_eff = min(rB, ya.to_pixel(0.02))
    return trace_curve_at(mask, xa, ya, (cL, rT, cR, rB_eff), at, half_w=3, inset_px=2, extra_err_px=1.2, search_w_fallback=10)


def read_fig5b_series(img: np.ndarray, conc_idx: int) -> list[dict]:
    """Read Fig. 5B bars across 8 construct conditions (conc_idx: 0=1uM white, 1=1.5uM gray, 2=2uM charcoal)."""
    H = img.shape[0]
    fr = find_plot_frame(img, row_range=(0, int(0.88 * H)))
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0, 5, 10, 15, 20, 25, 30, 35], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    g = gray(img)
    rB, rT = int(round(fr["rB"])), int(round(fr["rT"]))
    cL, cR = int(fr["left_band"][0]), int(round(fr["cR"])) - 4

    # Detect all 24 bars (8 groups x 3 bars: white, gray, dark charcoal)
    v_dark = (g[rB - 90 : rB - 10, cL:cR] < 130).mean(axis=0) > 0.70
    dark_runs = [(cL + s, cL + e) for s, e in runs_1d(v_dark)]

    bars = []  # (c0, c1, kind)
    thin = []
    for s, e in dark_runs:
        w = e - s + 1
        if 50 <= w <= 85:
            bars.append((s, e, "dark"))
        elif 2 <= w <= 14:
            thin.append(0.5 * (s + e))
    i = 0
    while i + 1 < len(thin):
        w = thin[i + 1] - thin[i]
        if 52 <= w <= 84:
            c0, c1 = int(round(thin[i])), int(round(thin[i + 1]))
            mg = float(np.median(g[rB - 80 : rB - 15, c0 + 10 : c1 - 10]))
            kind = "white" if mg > 215 else "gray"
            bars.append((c0, c1, kind))
            i += 2
        else:
            i += 1
    bars.sort(key=lambda b: b[0])
    if len(bars) != 24:
        raise ValueError(f"Fig. 5B: expected 24 bars, found {len(bars)}: {bars}")

    cal_err = max(ya.residual_px, 0.5)
    out = []
    for grp in range(8):
        c0, c1, kind = bars[3 * grp + conc_idx]
        c_mid = 0.5 * (c0 + c1) + 0.5
        q = list(range(c0 + 8, int(c_mid) - 10)) + list(range(int(c_mid) + 10, c1 - 7))
        if kind == "white":
            # Scan upward from bottom for the first horizontal dark border spanning >= 65% of q
            prof = (g[rT + 20 : rB - 30, q] < 135).mean(axis=1)
            hits = runs_1d(prof >= 0.65)
            if not hits:
                raise ValueError(f"Fig. 5B: no top border for white bar group {grp + 1}")
            s, e = hits[-1]
            r_top = (rT + 20) + 0.5 * (s + e) + 0.5
            fill_g = None
        else:
            r_cur = rB - 15
            while r_cur > rT + 20:
                if (g[r_cur, q] < 195).mean() < 0.5:
                    break
                r_cur -= 1
            r_top = float(r_cur + 1) + 0.5
            fill_g = float(np.median(g[int(r_top) + 12 : rB - 10, q]))
        y, lo, hi = val_interval(ya, r_top, 1.5 + cal_err)
        eb = measure_error_bar_vertical(
            img, c_mid, r_top, rT + 6, rB - 3, ya,
            dark_thresh=120.0, stem_half_w=4, gap_tol=3,
            mirror_lower_if_hidden=True, bar_fill_gray=fill_g,
        )
        out.append({
            "x": float(grp + 1), "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_mid, r_top], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"cond {grp + 1} ({kind}); {eb['note']}",
        })
    return out


def read_fig5c_bars(img: np.ndarray) -> list[dict]:
    """Read Fig. 5C (2 uM profilin-actin, 100 mM KCl, [VASP] = 0, 25, 50 nM) with SD error bars."""
    H = img.shape[0]
    fr = find_plot_frame(img, row_range=(0, int(0.78 * H)))
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0, 5, 10, 15, 20], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    g = gray(img)
    rB, rT = int(round(fr["rB"])), int(round(fr["rT"]))
    cL, cR = int(round(fr["cL"])) + 4, int(round(fr["cR"])) - 4

    # 3 bars: 1 white (border pair) + 2 gray (non-white runs)
    non_w = (g[rB - 90 : rB - 10, cL:cR] < 235).mean(axis=0) > 0.70
    gray_bars = [(cL + s, cL + e, "gray") for s, e in runs_1d(non_w) if e - s >= 80]
    # White bar is to the left of the first gray bar (two thin vertical dark borders ~120 px apart)
    v_d = (g[rB - 90 : rB - 10, cL : gray_bars[0][0] - 5] < 130).mean(axis=0) > 0.70
    thin = [cL + 0.5 * (s + e) for s, e in runs_1d(v_d) if e - s <= 12]
    if len(thin) < 2 or len(gray_bars) != 2:
        raise ValueError(f"Fig. 5C: unexpected bars thin={thin} gray={gray_bars}")
    bars = [(int(round(thin[0])), int(round(thin[-1])), "white")] + gray_bars
    concs = [0.0, 25.0, 50.0]
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for conc, (c0, c1, kind) in zip(concs, bars):
        c_mid = 0.5 * (c0 + c1) + 0.5
        q = list(range(c0 + 12, int(c_mid) - 14)) + list(range(int(c_mid) + 14, c1 - 12))
        if kind == "white":
            prof = (g[rT + 30 : rB - 30, q] < 135).mean(axis=1)
            s, e = runs_1d(prof >= 0.65)[-1]
            r_top = (rT + 30) + 0.5 * (s + e) + 0.5
            fill_g = None
        else:
            r_cur = rB - 20
            while r_cur > rT + 10:
                if (g[r_cur, q] < 195).mean() < 0.5:
                    break
                r_cur -= 1
            r_top = float(r_cur + 1) + 0.5
            fill_g = float(np.median(g[int(r_top) + 15 : rB - 15, q]))
        y, lo, hi = val_interval(ya, r_top, 1.5 + cal_err)
        eb = measure_error_bar_vertical(
            img, c_mid, r_top, rT + 4, rB - 4, ya,
            dark_thresh=115.0, stem_half_w=4, gap_tol=3,
            mirror_lower_if_hidden=True, bar_fill_gray=fill_g,
        )
        out.append({
            "x": conc, "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_mid, r_top], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"Fig. 5C {kind}; {eb['note']}",
        })
    return out


def read_fig5de_bars(img: np.ndarray, y_ticks_vals: list[float], group: str) -> list[dict]:
    """Read Fig. 5D or Fig. 5E: 4 gray bars -> hProI (bars 0,1 at 1,2 uM) or hProI_H133S (bars 2,3 at 1,2 uM)."""
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, y_ticks_vals, increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    g = gray(img)
    rB, rT = int(round(fr["rB"])), int(round(fr["rT"]))
    cL, cR = int(round(fr["cL"])) + 4, int(round(fr["cR"])) - 4
    non_w = (g[rB - 90 : rB - 10, cL:cR] < 235).mean(axis=0) > 0.70
    bcols = [(cL + s, cL + e) for s, e in runs_1d(non_w) if e - s >= 80]
    if len(bcols) != 4:
        raise ValueError(f"Fig. 5D/E: expected 4 bars, found {len(bcols)}: {bcols}")
    idxs = [0, 1] if group == "hpro1" else [2, 3]
    concs = [1.0, 2.0]
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for conc, bi in zip(concs, idxs):
        c0, c1 = bcols[bi]
        c_mid = 0.5 * (c0 + c1) + 0.5
        q = list(range(c0 + 12, int(c_mid) - 14)) + list(range(int(c_mid) + 14, c1 - 12))
        r_cur = rB - 20
        while r_cur > rT + 10:
            if (g[r_cur, q] < 210).mean() < 0.5:
                break
            r_cur -= 1
        r_top = float(r_cur + 1) + 0.5
        fill_g = float(np.median(g[int(r_top) + 12 : rB - 12, q]))
        y, lo, hi = val_interval(ya, r_top, 1.5 + cal_err)
        eb = measure_error_bar_vertical(
            img, c_mid, r_top, rT + 4, rB - 4, ya,
            dark_thresh=115.0, stem_half_w=4, gap_tol=3,
            mirror_lower_if_hidden=True, bar_fill_gray=fill_g,
        )
        out.append({
            "x": conc, "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_mid, r_top], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"{group} {conc:g} uM; {eb['note']}",
        })
    return out


# ==========================================================================
# Batch 4 helpers: Figs. 6B, 6C, 6D, S2D(i-iiii)
# ==========================================================================


def _x_axis_from_labels(img: np.ndarray, fr: dict, values: list[float], gap_merge: int = 18) -> D.Axis:
    g = gray(img)
    rB = int(round(fr["x_row"]))
    c0, c1 = max(0, int(fr["y_col"]) - 35), min(img.shape[1], int(fr["right_col"]) + 35)
    sub = g[rB + 10 : min(img.shape[0], rB + 58), c0:c1] < 160
    col_has = ndimage.binary_closing(sub.any(axis=0), structure=np.ones(gap_merge, bool))
    lbls = [c0 + 0.5 * (s + e) + 0.5 for s, e in runs_1d(col_has) if e - s >= 4]
    return match_linear_axis(lbls, values, increasing_px=True, tol_px=14.0, max_resid=8.0)


def calibrate_fig6b(img: np.ndarray):
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000], increasing_px=False)
    xa = _x_axis_from_labels(img, fr, [0, 5, 10, 15, 20], gap_merge=18)
    return xa, ya, fr


def calibrate_fig6c(img: np.ndarray):
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False, extra_anchors=[fr["x_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 600.0)])
    return xa, ya, fr


def calibrate_fig6d(img: np.ndarray):
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 25.0)])
    return xa, ya, fr


def calibrate_figS2d_i(img: np.ndarray):
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0, 1000, 2000, 3000, 4000, 5000, 6000], increasing_px=False, extra_anchors=[fr["top_row"]])
    xa = _x_axis_from_labels(img, fr, [0, 2, 4, 6, 8, 10, 12], gap_merge=14)
    return xa, ya, fr


def calibrate_figS2d_ii(img: np.ndarray):
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 10.0)])
    return xa, ya, fr


def calibrate_figS2d_iii(img: np.ndarray):
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", thresh=190.0, probe=(2, 10), spine_band=fr["left_band"])
    ya = match_linear_axis(yi, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], increasing_px=False, extra_anchors=[fr["x_row"]])
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 500.0)])
    return xa, ya, fr


def read_ranked_dark_curves(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    fr: dict,
    at: list[float],
    rank_fn: Callable[[float, list[tuple[float, float, float]]], int | None],
    mask: np.ndarray | None = None,
    dark_thresh: float = 135.0,
    half_w: int = 2,
    merge_gap: int = 6,
    y_min_val: float = -0.02,
    y_max_val: float = 1.04,
) -> list[dict]:
    """Read one curve out of a vertically ordered family of curves (e.g. Fig. 6B/C, S2D, S5C).

    `rank_fn(x, clusters)` receives the top-to-bottom sorted list of `(r_center, r_top, r_bot)`
    clusters in the column window and returns the index to pick (or None if missing).
    """
    m = (gray(img) < dark_thresh) if mask is None else mask
    cL, rT, cR, rB = fr["box"]
    r_lo = max(int(math.ceil(rT + 2)), int(math.floor(ya.to_pixel(y_max_val))))
    r_hi = min(int(math.floor(rB - 2)), int(math.ceil(ya.to_pixel(y_min_val))))
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for x in at:
        xv = float(x)
        c_f = xa.to_pixel(xv)
        c0 = max(int(math.ceil(cL + 1)), int(round(c_f)) - half_w)
        c1 = min(int(math.floor(cR - 1)), int(round(c_f)) + half_w + 1)
        sub = m[r_lo:r_hi, c0:c1]
        row_any = sub.any(axis=1)
        if merge_gap > 1:
            row_any = ndimage.binary_closing(row_any, structure=np.ones(merge_gap, bool))
        clusters = []
        for s, e in runs_1d(row_any):
            ys = np.nonzero(sub[s : e + 1])[0]
            if len(ys) == 0:
                continue
            rc = r_lo + s + float(np.median(ys)) + 0.5
            clusters.append((rc, float(r_lo + s), float(r_lo + e + 1)))
        idx = rank_fn(xv, clusters)
        if idx is None or not (-len(clusters) <= idx < len(clusters)):
            out.append({"x": xv, "y": None, "note": f"expected curve cluster not found (n={len(clusters)})"})
            continue
        rc, rt, rb = clusters[idx]
        half_px = 0.5 * max(rb - rt, 2.0) + cal_err + 1.0
        y, lo, hi = val_interval(ya, rc, half_px)
        out.append({
            "x": xv, "y": y, "lo": lo, "hi": hi,
            "px": [float(c_f), float(rc)],
            "px_lo": float(rc - half_px), "px_hi": float(rc + half_px),
            "note": f"cluster {idx}/{len(clusters)} (thick {rb - rt:.1f}px)",
        })
    return out


# ==========================================================================
# Batch 5 helpers: Figs. S1B, S3D, S3E, S4D, S5C
# ==========================================================================


def read_figS1b_bars(img: np.ndarray) -> list[dict]:
    """Read Fig. S1B horizontal bar chart (5 conditions at 1 uM actin, SD error bars)."""
    fr = find_plot_frame(img)
    xt = detect_spine_ticks(img, fr["top_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "out", spine_band=fr["top_band"])
    # Note: top axis measures polymerization rate (0..35 sub/sec)
    rate_ax = match_linear_axis(xt, [0, 5, 10, 15, 20, 25, 30, 35], increasing_px=True, extra_anchors=[fr["y_col"], fr["right_col"]])
    g = gray(img)
    cL, cR = int(round(fr["cL"])) + 4, int(round(fr["cR"])) - 4
    rT, rB = int(round(fr["rT"])) + 4, int(round(fr["rB"])) - 4
    nonw_h = (g[rT:rB, cL + 10 : cL + 120] < 220).mean(axis=1) > 0.70
    hbars = [(rT + s, rT + e) for s, e in runs_1d(nonw_h) if e - s >= 50]
    if len(hbars) != 5:
        raise ValueError(f"Fig. S1B: expected 5 horizontal bars, found {len(hbars)}: {hbars}")
    cal_err = max(rate_ax.residual_px, 0.5)
    labels = ["1uM Actin", "+ VASP", "+ VASP (LIL-AAA)", "+ VASP (RK-EE)", "+ VASP (RRRK-4E)"]
    out = []
    for idx, ((r0, r1), lbl) in enumerate(zip(hbars, labels), start=1):
        r_mid = 0.5 * (r0 + r1) + 0.5
        # Use upper and lower quarters of bar height to avoid horizontal center error-bar stem
        q_rows = list(range(r0 + 10, int(r_mid) - 12)) + list(range(int(r_mid) + 12, r1 - 10))
        # Scan rightward from cL + 15 until gray bar fill ends (turns white > 225)
        c_cur = cL + 15
        while c_cur < cR - 5:
            if (g[q_rows, c_cur] < 210).mean() < 0.5:
                break
            c_cur += 1
        c_end = float(c_cur - 1) + 0.5
        fill_g = float(np.median(g[q_rows, cL + 15 : int(c_end) - 10]))
        y, lo, hi = val_interval(rate_ax, c_end, 1.5 + cal_err)
        # Trace horizontal error bar rightward and leftward from c_end along r_mid +- 4
        rm0, rm1 = int(round(r_mid)) - 4, int(round(r_mid)) + 5
        # Right (upper) error bar
        c_hi = int(math.ceil(c_end))
        best_r, gap = c_hi, 0
        for c in range(c_hi + 1, cR):
            if (g[rm0:rm1, c] < 115.0).any():
                best_r, gap = c, 0
            else:
                gap += 1
                if gap > 3:
                    break
        # Left (lower) error bar inside gray bar (stricter threshold darker than gray fill)
        in_thr = min(105.0, fill_g - 35.0)
        best_l, gap = int(math.floor(c_end)), 0
        for c in range(int(math.floor(c_end)) - 4, cL + 5, -1):
            if (g[rm0:rm1, c] < in_thr).any():
                best_l, gap = c, 0
            else:
                gap += 1
                if gap > 3:
                    break
        err_px = 1.5 + cal_err
        bar_hi = float(rate_ax.to_value(best_r + 0.5))
        bar_hi_out = float(rate_ax.to_value(best_r + 0.5 + err_px))
        if best_l < int(math.floor(c_end)) - 5:
            bar_lo = float(rate_ax.to_value(best_l + 0.5))
            bar_lo_out = float(rate_ax.to_value(best_l + 0.5 - err_px))
            enote = f"horiz err bar [{best_l}, {best_r}]px"
        else:
            bar_lo = 2.0 * y - bar_hi
            bar_lo_out = 2.0 * y - bar_hi_out
            enote = f"right err bar {best_r}px (left mirrored)"
        out.append({
            "x": float(idx), "y": y, "lo": lo, "hi": hi,
            "bar_lo": min(bar_lo, bar_hi), "bar_hi": max(bar_lo, bar_hi),
            "bar_lo_outer": min(bar_lo_out, bar_hi_out), "bar_hi_outer": max(bar_lo_out, bar_hi_out),
            "px": [c_end, r_mid],
            "note": f"{lbl}; {enote}",
        })
    return out


def read_vertical_gray_bars(
    img: np.ndarray,
    y_ticks_vals: list[float],
    n_expected: int,
    labels: list[str],
    min_bar_w: int = 60,
    extra_y_anchors: bool = True,
) -> list[dict]:
    """Generic calibrated reader for clean vertical gray bar charts with black error bars (Figs. S3D, S4D)."""
    fr = find_plot_frame(img)
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    anchors = [fr["x_row"], fr["top_row"]] if extra_y_anchors else [fr["x_row"]]
    ya = match_linear_axis(yi, y_ticks_vals, increasing_px=False, extra_anchors=anchors)
    g = gray(img)
    rB, rT = int(round(fr["rB"])), int(round(fr["rT"]))
    cL, cR = int(round(fr["cL"])) + 4, int(round(fr["cR"])) - 4
    nonw = (g[rB - 90 : rB - 10, cL:cR] < 225).mean(axis=0) > 0.70
    bcols = [(cL + s, cL + e) for s, e in runs_1d(nonw) if e - s >= min_bar_w]
    if len(bcols) != n_expected:
        raise ValueError(f"Expected {n_expected} gray bars, found {len(bcols)}: {bcols}")
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for idx, ((c0, c1), lbl) in enumerate(zip(bcols, labels), start=1):
        c_mid = 0.5 * (c0 + c1) + 0.5
        q = list(range(c0 + 14, int(c_mid) - 16)) + list(range(int(c_mid) + 16, c1 - 14))
        r_cur = rB - 15
        while r_cur > rT + 10:
            if (g[r_cur, q] < 205).mean() < 0.5:
                break
            r_cur -= 1
        r_top = float(r_cur + 1) + 0.5
        fill_g = float(np.median(g[int(r_top) + 14 : rB - 12, q]))
        y, lo, hi = val_interval(ya, r_top, 1.5 + cal_err)
        eb = measure_error_bar_vertical(
            img, c_mid, r_top, rT + 4, rB - 4, ya,
            dark_thresh=110.0, stem_half_w=4, gap_tol=3,
            mirror_lower_if_hidden=True, bar_fill_gray=fill_g,
        )
        out.append({
            "x": float(idx), "y": y, "lo": lo, "hi": hi,
            "bar_lo": eb["bar_lo"], "bar_hi": eb["bar_hi"],
            "bar_lo_outer": eb["bar_lo_outer"], "bar_hi_outer": eb["bar_hi_outer"],
            "px": [c_mid, r_top], "px_lo": eb["px_lo"], "px_hi": eb["px_hi"],
            "note": f"{lbl}; {eb['note']}",
        })
    return out


def calibrate_figS3e(img: np.ndarray):
    fr = find_plot_frame(img)
    xa = D.Axis("linear", [(fr["y_col"], 0.0), (fr["right_col"], 2400.0)])
    ya = D.Axis("linear", [(fr["x_row"], 0.0), (fr["top_row"], 1.0)])
    return xa, ya, fr


def calibrate_figS5c(img: np.ndarray):
    fr = find_plot_frame(img, min_h_frac=0.30, min_v_frac=0.30)
    xi = detect_spine_ticks(img, fr["x_row"], int(fr["cL"]) - 5, int(fr["cR"]) + 5, "x", "in", spine_band=fr["bot_band"])
    yi = detect_spine_ticks(img, fr["y_col"], int(fr["rT"]) - 5, int(fr["rB"]) + 5, "y", "in", spine_band=fr["left_band"])
    xa = match_linear_axis(xi, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8], increasing_px=True, extra_anchors=[fr["y_col"], fr["right_col"]])
    ya = match_linear_axis(yi, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0], increasing_px=False, extra_anchors=[fr["x_row"], fr["top_row"]])
    return xa, ya, fr


def detect_open_circles(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    fr: dict,
    exclude_fn: Callable[[float, float], bool] | None = None,
    thr: float = 232.0,
    close_k: int = 7,
    area_range: tuple[int, int] = (80, 620),
    dim_range: tuple[int, int] = (10, 34),
) -> list[dict]:
    """Detect hollow open-circle plot markers via ring closing + enclosed white-hole centroids."""
    g = gray(img)
    cL, rT, cR, rB = [int(round(v)) for v in fr["box"]]
    ring = ndimage.binary_closing(g < thr, structure=np.ones((close_k, close_k), bool))
    holes = ~ring
    holes[: rT + 2, :] = False
    holes[rB - 2 :, :] = False
    holes[:, : cL + 2] = False
    holes[:, cR - 2 :] = False
    lab, n = ndimage.label(holes)
    # Remove any component touching the box interior boundary (background white)
    border_ids = (
        set(np.unique(lab[rT + 2, cL + 2 : cR - 2]))
        | set(np.unique(lab[rB - 3, cL + 2 : cR - 2]))
        | set(np.unique(lab[rT + 2 : rB - 2, cL + 2]))
        | set(np.unique(lab[rT + 2 : rB - 2, cR - 3]))
    )
    cal_err = max(ya.residual_px, 0.5)
    pts = []
    for i in range(1, n + 1):
        if i in border_ids:
            continue
        rr, cc = np.nonzero(lab == i)
        area = len(rr)
        if not (area_range[0] <= area <= area_range[1]):
            continue
        h = int(rr.max() - rr.min() + 1)
        w = int(cc.max() - cc.min() + 1)
        if not (dim_range[0] <= h <= dim_range[1] and dim_range[0] <= w <= dim_range[1]):
            continue
        if not (0.62 <= h / w <= 1.62):
            continue
        cm, rm = float(cc.mean()) + 0.5, float(rr.mean()) + 0.5
        xv = float(xa.to_value(cm))
        if exclude_fn is not None and exclude_fn(xv, float(ya.to_value(rm))):
            continue
        y, lo, hi = val_interval(ya, rm, 2.0 + cal_err)
        pts.append({
            "x": round(xv, 3), "y": y, "lo": lo, "hi": hi,
            "px": [cm, rm], "px_lo": rm - 2.0 - cal_err, "px_hi": rm + 2.0 + cal_err,
            "note": f"open circle hole {area}px ({w}x{h})",
        })
    pts.sort(key=lambda p: p["x"])
    return pts


def trace_monotonic_capping_curve(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    fr: dict,
    seed_x: float,
    seed_rank: int,
    at: list[float],
    mask_fn: Callable[[np.ndarray, D.Axis, D.Axis], None] | None = None,
    seed_thresh: float = 70.0,
    march_thresh: float = 95.0,
    max_gap: int = 60,
    y_min: float = 0.02,
    y_max: float = 1.02,
    strict_mono_tol: float = 1.0,
) -> list[dict]:
    """Trace one smooth, monotonically increasing black fit curve from a clean seed column.

    Enforces `d(row)/d(col) <= 0` (monotonically non-decreasing y) in both left and right
    marching directions so downward-right callout leader lines can never capture the tracker.
    """
    g = gray(img)
    m_seed = g < seed_thresh
    m_march = g < march_thresh
    if mask_fn is not None:
        mask_fn(m_seed, xa, ya)
        mask_fn(m_march, xa, ya)

    cL, rT, cR, rB = fr["box"]
    r_lo = max(int(math.ceil(rT + 2)), int(math.floor(ya.to_pixel(y_max))))
    r_hi = min(int(math.floor(rB - 2)), int(math.ceil(ya.to_pixel(y_min))))
    c_min = int(math.ceil(cL + 2))
    c_max = int(math.floor(cR - 2))
    cs = int(round(xa.to_pixel(seed_x)))

    sub = m_seed[r_lo:r_hi, cs - 1 : cs + 2]
    row_any = ndimage.binary_closing(sub.any(axis=1), structure=np.ones(5, bool))
    seed_cls = []
    for s, e in runs_1d(row_any):
        ys = np.nonzero(sub[s : e + 1])[0]
        if len(ys):
            seed_cls.append(r_lo + s + float(np.median(ys)) + 0.5)
    if not (-len(seed_cls) <= seed_rank < len(seed_cls)):
        raise ValueError(f"seed x={seed_x}: rank {seed_rank} out of {len(seed_cls)} clusters")
    r_seed = seed_cls[seed_rank]

    def march(direction: int) -> dict[int, tuple[float, float]]:
        pos: dict[int, tuple[float, float]] = {cs: (r_seed, 4.0)}
        hist: list[tuple[int, float]] = [(cs, r_seed)]
        c = cs
        gap = 0
        while True:
            c += direction
            if c < c_min or c > c_max:
                break
            recent = hist[-30:]
            if len(recent) >= 6 and recent[-1][0] != recent[0][0]:
                slope = (recent[-1][1] - recent[0][1]) / (recent[-1][0] - recent[0][0])
                slope = min(0.0, max(-12.0, slope))
            else:
                slope = -0.35
            r_prev = hist[-1][1]
            dc = c - hist[-1][0]
            r_pred = r_prev + slope * dc
            max_step = max(14.0, abs(slope * dc) * 2.0 + 8.0)
            if direction > 0:
                w_lo = max(r_lo, int(math.floor(r_pred - max_step)))
                w_hi = min(r_hi, int(math.ceil(r_prev + max(strict_mono_tol, 1.5))))
            else:
                w_lo = max(r_lo, int(math.floor(r_prev - strict_mono_tol)))
                w_hi = min(r_hi, int(math.ceil(r_pred + max_step)))
            if w_lo >= w_hi:
                gap += 1
                if gap > max_gap:
                    break
                continue
            runs = runs_1d(ndimage.binary_closing(m_march[w_lo:w_hi, c], structure=np.ones(3, bool)))
            best = None
            for s, e in runs:
                thick = float(e - s + 1)
                if thick > max(55.0, abs(slope) * 12.0 + 22.0):
                    continue
                rc = w_lo + 0.5 * (s + e) + 0.5
                dist = abs(rc - r_pred)
                if best is None or dist < best[0]:
                    best = (dist, rc, thick)
            if best is not None and best[0] <= max(14.0, abs(slope * dc) * 1.4 + 10.0):
                pos[c] = (best[1], best[2])
                hist.append((c, best[1]))
                gap = 0
            else:
                gap += 1
                if gap > max_gap:
                    break
        return pos

    pos = {cs: (r_seed, 4.0)}
    pos.update(march(+1))
    pos.update(march(-1))

    # Interpolate across tiny bridged gaps so any query column between tracked endpoints is covered
    tracked_cols = sorted(pos.keys())
    cal_err = max(ya.residual_px, 0.5)
    out = []
    for x in at:
        xv = float(x)
        c_f = xa.to_pixel(xv)
        ci = int(round(c_f))
        near = [pos[cc] for cc in range(ci - 3, ci + 4) if cc in pos]
        if near:
            rc = float(np.median([p[0] for p in near]))
            thick = float(np.median([p[1] for p in near]))
            extra = 1.0
            note_tag = f"traced (thick {thick:.1f}px)"
        else:
            # Linear interpolation across a short bridged gap (< max_gap px)
            lefts = [cc for cc in tracked_cols if cc < ci]
            rights = [cc for cc in tracked_cols if cc > ci]
            if lefts and rights and (rights[0] - lefts[-1]) <= max_gap + 5:
                cl, cr = lefts[-1], rights[0]
                t_frac = (ci - cl) / (cr - cl)
                rc = pos[cl][0] * (1.0 - t_frac) + pos[cr][0] * t_frac
                thick = max(pos[cl][1], pos[cr][1])
                extra = 2.5
                note_tag = f"interpolated across {cr - cl}px gap"
            else:
                out.append({"x": xv, "y": None, "note": f"traced curve did not reach x={xv}"})
                continue
        half_px = 0.5 * max(thick, 3.0) + cal_err + extra
        y, lo, hi = val_interval(ya, rc, half_px)
        out.append({
            "x": xv, "y": y, "lo": lo, "hi": hi,
            "px": [float(c_f), float(rc)],
            "px_lo": float(rc - half_px), "px_hi": float(rc + half_px),
            "note": note_tag,
        })
    return out


def read_figS5c_speed(img: np.ndarray, at: list[float], speed_rpm: int, x_cross: float) -> list[dict]:
    """Read one of the 3 red best-fit AUC curves (5000, 7000, 10000 RPM) in Fig. S5C.

    Left of the mutual crossing `x_cross`, vertical order top->bottom is [5000, 7000, 10000].
    Right of `x_cross`, vertical order top->bottom reverses to [10000, 7000, 5000].
    """
    xa, ya, fr = calibrate_figS5c(img)
    rm = red_mask(img, r_min=115.0, rg_diff=35.0, rb_diff=35.0)
    rank_left = {5000: 0, 7000: 1, 10000: 2}[speed_rpm]
    rank_right = {10000: 0, 7000: 1, 5000: 2}[speed_rpm]

    def pick(xv, clusters):
        # Filter out tiny < 4 px stray red artifacts
        valid = [i for i, (_, rt, rb) in enumerate(clusters) if (rb - rt) >= 4.0]
        if len(valid) == 3:
            k = rank_left if xv < x_cross else rank_right
            return valid[k]
        return None

    return read_ranked_dark_curves(
        img, xa, ya, fr, at, pick,
        mask=rm, half_w=1, merge_gap=5, y_min_val=-0.1, y_max_val=12.1,
    )
