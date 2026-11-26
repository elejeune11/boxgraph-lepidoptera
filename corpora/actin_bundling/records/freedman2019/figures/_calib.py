"""Shared image-based axis calibration and data extraction helpers for freedman2019 figures.

All functions operate purely on the passed RGB `image` numpy array (no PDF rendering or
external I/O). Axis frames, tick marks, markers, error bars, bars, and curves are detected
directly from the image pixels.
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np
from scipy import ndimage

from boxgraph import digitize as D


# ---------------------------------------------------------------------------
# Basic color & geometry helpers
# ---------------------------------------------------------------------------

def _rgb(img: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    f = img.astype(float)
    return f[:, :, 0], f[:, :, 1], f[:, :, 2]


def black_mask(img: np.ndarray, thresh: float = 120.0, sat_max: float = 40.0) -> np.ndarray:
    r, g, b = _rgb(img)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    return (r < thresh) & (g < thresh) & (b < thresh) & (sat < sat_max)


def red_mask(img: np.ndarray, r_min: float = 145.0, diff: float = 45.0) -> np.ndarray:
    r, g, b = _rgb(img)
    return (r > r_min) & (r - g > diff) & (r - b > diff)


def cyan_mask(img: np.ndarray, gb_min: float = 130.0, diff: float = 40.0) -> np.ndarray:
    r, g, b = _rgb(img)
    return (g > gb_min) & (b > gb_min) & (g - r > diff) & (b - r > diff)


def rgb_dist_mask(img: np.ndarray, target: Sequence[float], tol: float = 40.0) -> np.ndarray:
    f = img.astype(float)
    t = np.asarray(target, dtype=float)
    return np.sqrt(((f - t) ** 2).sum(axis=2)) <= tol


def _runs_1d(arr: np.ndarray, max_gap: int = 1) -> list[tuple[int, int]]:
    idx = np.nonzero(arr)[0]
    if len(idx) == 0:
        return []
    br = np.where(np.diff(idx) > max_gap)[0]
    starts = np.r_[idx[0], idx[br + 1]]
    ends = np.r_[idx[br], idx[-1]]
    return [(int(s), int(e)) for s, e in zip(starts, ends)]


def find_frame_spines(img: np.ndarray, min_frac: float = 0.35) -> tuple[list[tuple[int, int]], list[tuple[int, int]], np.ndarray]:
    """Locate closed black frame horizontal and vertical spine bands [(start, end), ...]."""
    blk = black_mask(img, thresh=120.0, sat_max=40.0)
    H, W = blk.shape

    def longest_run(a: np.ndarray) -> int:
        rs = _runs_1d(a, max_gap=2)
        return max((e - s + 1 for s, e in rs), default=0)

    hs = np.array([r for r in range(H) if blk[r].sum() > min_frac * W and longest_run(blk[r]) > min_frac * W])
    vs = np.array([c for c in range(W) if blk[:, c].sum() > min_frac * H and longest_run(blk[:, c]) > min_frac * H])
    if len(hs) < 2 or len(vs) < 2:
        raise ValueError(f"Could not locate closed plot frame (H={len(hs)}, V={len(vs)})")

    def bands(idx: np.ndarray) -> list[tuple[int, int]]:
        br = np.where(np.diff(idx) > 3)[0]
        s = np.r_[idx[0], idx[br + 1]]
        e = np.r_[idx[br], idx[-1]]
        return [(int(a), int(b)) for a, b in zip(s, e)]

    return bands(hs), bands(vs), blk


def _detect_ticks_strip(strip: np.ndarray, min_len: int = 5, max_width: int = 20) -> list[tuple[float, int]]:
    """Detect ticks protruding from row 0 of a boolean strip (depth x along_axis).

    Returns list of (center_offset, max_depth).
    """
    attached = strip[:3].mean(axis=0) > 0.5
    depths = np.zeros(strip.shape[1], dtype=int)
    for j in np.nonzero(attached)[0]:
        col = strip[:, j]
        k = 0
        while k < len(col) and col[k]:
            k += 1
        depths[j] = k
    good = depths >= min_len
    out: list[tuple[float, int]] = []
    for s, e in _runs_1d(good, max_gap=2):
        if e - s + 1 <= max_width:
            out.append((0.5 * (s + e), int(depths[s : e + 1].max())))
    return out


def find_inward_ticks(
    blk: np.ndarray,
    rT: tuple[int, int],
    rB: tuple[int, int],
    cL: tuple[int, int],
    cR: tuple[int, int],
    side: str,
    depth: int = 45,
    min_len: int = 5,
    max_depth: int = 35,
) -> list[tuple[float, int]]:
    """Find inward ticks on Mathematica-style closed frame."""
    if side == "bottom":
        strip = blk[rB[0] - 1 : max(0, rB[0] - 1 - depth) : -1, cL[1] + 2 : cR[0] - 2]
        offset = cL[1] + 2
    elif side == "top":
        strip = blk[rT[1] + 1 : min(blk.shape[0], rT[1] + 1 + depth), cL[1] + 2 : cR[0] - 2]
        offset = cL[1] + 2
    elif side == "left":
        strip = blk[rT[1] + 2 : rB[0] - 2, cL[1] + 1 : min(blk.shape[1], cL[1] + 1 + depth)].T
        offset = rT[1] + 2
    elif side == "right":
        strip = blk[rT[1] + 2 : rB[0] - 2, cR[0] - 1 : max(0, cR[0] - 1 - depth) : -1].T
        offset = rT[1] + 2
    else:
        raise ValueError(f"Unknown side {side!r}")
    raw = _detect_ticks_strip(strip, min_len=min_len)
    return [(offset + p, d) for p, d in raw if d <= max_depth]


def find_outward_ticks(
    blk: np.ndarray,
    rT: tuple[int, int],
    rB: tuple[int, int],
    cL: tuple[int, int],
    cR: tuple[int, int],
    side: str,
    depth: int = 65,
    min_len: int = 8,
    max_depth: int = 58,
) -> list[tuple[float, int]]:
    """Find outward ticks on Matplotlib-style closed frame."""
    H, W = blk.shape
    if side == "bottom":
        strip = blk[rB[1] + 1 : min(H, rB[1] + 1 + depth), cL[0] : cR[1] + 1]
        offset = cL[0]
    elif side == "left":
        strip = blk[rT[0] : rB[1] + 1, max(0, cL[0] - depth) : cL[0]][:, ::-1].T
        offset = rT[0]
    elif side == "right":
        strip = blk[rT[0] : rB[1] + 1, cR[1] + 1 : min(W, cR[1] + 1 + depth)].T
        offset = rT[0]
    else:
        raise ValueError(f"Unknown side {side!r}")
    raw = _detect_ticks_strip(strip, min_len=min_len)
    out = []
    for p, d in raw:
        pos = offset + p
        if d <= max_depth:
            out.append((pos, d))
        elif side in ("left", "right") and abs(pos - 0.5 * (rB[0] + rB[1])) <= 4.0:
            # Shared bottom spine extending into adjacent subplot also sits at y=0 tick
            out.append((pos, 45))
    return out


def pick_major_ticks(ticks: list[tuple[float, int]], expected_n: int, min_depth_frac: float = 0.72) -> list[float]:
    """Select major ticks by depth threshold, falling back to longest expected_n ticks if needed."""
    if not ticks:
        raise ValueError("No ticks detected")
    depths = sorted(d for _, d in ticks)
    ref_d = depths[min(len(depths) - 1, max(0, len(depths) - expected_n))]
    maj = [p for p, d in ticks if d >= ref_d * min_depth_frac]
    if len(maj) == expected_n:
        return sorted(maj)
    uniq_d = sorted(set(depths), reverse=True)
    for d_cut in uniq_d:
        cand = sorted(p for p, d in ticks if d >= d_cut - 2 and d <= d_cut + 3)
        if len(cand) == expected_n:
            return cand
    top = sorted(ticks, key=lambda t: t[1], reverse=True)[:expected_n]
    return sorted(p for p, _ in top)


def match_ticks_regular(ticks: list[tuple[float, int]], values: Sequence[float], scale: str = "linear") -> D.Axis:
    """Calibrate an axis from detected ticks by matching the expected tick values.

    Uses major ticks first; verifies low residual (< 2.0 px).
    """
    n = len(values)
    pos = pick_major_ticks(ticks, n)
    ax = D.Axis(scale, list(zip(pos, values)))
    if ax.residual_px <= 2.0:
        return ax
    all_p = sorted(p for p, _ in ticks)
    if len(all_p) >= n:
        best_ax = ax
        u_vals = np.array([math.log10(v) if scale == "log10" else float(v) for v in values])
        u_norm = (u_vals - u_vals[0]) / (u_vals[-1] - u_vals[0])
        for i in range(len(all_p)):
            for j in range(i + n - 1, len(all_p)):
                p0, p1 = all_p[i], all_p[j]
                pred = p0 + u_norm * (p1 - p0)
                matched = []
                ok = True
                for pr in pred:
                    nearest = min(all_p, key=lambda x: abs(x - pr))
                    if abs(nearest - pr) > 4.0:
                        ok = False
                        break
                    matched.append(nearest)
                if ok and len(set(matched)) == n:
                    cand_ax = D.Axis(scale, list(zip(matched, values)))
                    if cand_ax.residual_px < best_ax.residual_px:
                        best_ax = cand_ax
        if best_ax.residual_px <= 2.0:
            return best_ax
    raise ValueError(f"Tick calibration misfit too large ({ax.residual_px:.2f} px) for values {values}")


def _make_row(
    x: float,
    y: float,
    col: float,
    row: float,
    ya: D.Axis,
    half_px: float,
    bar_lo_row: float | None = None,
    bar_hi_row: float | None = None,
    bar_pad_px: float = 2.0,
    note: str = "",
) -> dict:
    """Construct a BoxGraph code-mode reading row with calibrated interval and optional error bars."""
    y_hi = ya.to_value(row - half_px)
    y_lo = ya.to_value(row + half_px)
    lo, hi = min(y_lo, y_hi), max(y_lo, y_hi)
    lo = min(lo, y)
    hi = max(hi, y)
    rec: dict = {
        "x": float(x),
        "y": float(y),
        "lo": float(lo),
        "hi": float(hi),
        "px": [float(col), float(row)],
        "px_lo": float(row - half_px),
        "px_hi": float(row + half_px),
        "note": note,
    }
    if bar_lo_row is not None and bar_hi_row is not None:
        top_r = min(bar_hi_row, row)
        bot_r = max(bar_lo_row, row)
        b_hi = ya.to_value(top_r)
        b_lo = ya.to_value(bot_r)
        b_hi_out = ya.to_value(top_r - bar_pad_px)
        b_lo_out = ya.to_value(bot_r + bar_pad_px)
        b_lo, b_hi = min(b_lo, y), max(b_hi, y)
        b_lo_out, b_hi_out = min(b_lo_out, lo, b_lo), max(b_hi_out, hi, b_hi)
        rec["bar_lo"] = float(b_lo)
        rec["bar_hi"] = float(b_hi)
        rec["bar_lo_outer"] = float(b_lo_out)
        rec["bar_hi_outer"] = float(b_hi_out)
        rec["px_lo"] = float(top_r)
        rec["px_hi"] = float(bot_r)
    return rec


# ---------------------------------------------------------------------------
# Archetype 1: Mathematica marker + error-bar plots (inward ticks)
# ---------------------------------------------------------------------------

def _measure_cap_and_center(
    mask: np.ndarray,
    col_c: int,
    r0: int,
    r1: int,
    cap_half_w: int = 34,
    cap_frac: float = 0.65,
    y_hint_range: tuple[int, int] | None = None,
) -> tuple[float, float, float, str]:
    """At a marker column `col_c`, locate horizontal error-bar caps via bilateral ink-count detection.

    A Mathematica error-bar cap spans roughly +/- 32-34 px about `col_c` (total width ~64-69 px).
    Requiring substantial ink on BOTH left (`col_c - 34 .. col_c - 4`) and right
    (`col_c + 4 .. col_c + 34`) sides rejects vertical stems (~10 px wide), one-sided connecting
    lines, and thin open-marker walls while staying robust even when an overlapping marker or
    dashed line slices part of a cap.
    Returns (center_row, top_cap_row, bot_cap_row, note).
    """
    H, W = mask.shape
    rr0 = max(r0, y_hint_range[0]) if y_hint_range else r0
    rr1 = min(r1, y_hint_range[1]) if y_hint_range else r1

    hw = max(32, cap_half_w)
    lc0, lc1 = max(0, col_c - hw), max(0, col_c - 4)
    rc0, rc1 = min(W, col_c + 5), min(W, col_c + hw + 1)
    fc0, fc1 = max(0, col_c - hw), min(W, col_c + hw + 1)

    left_cnt = mask[rr0:rr1, lc0:lc1].sum(axis=1)
    right_cnt = mask[rr0:rr1, rc0:rc1].sum(axis=1)
    tot_cnt = mask[rr0:rr1, fc0:fc1].sum(axis=1)

    # First try strict cap criterion (wide bilateral cap, >= 42 px total and >= 12 px each side)
    # Then relax slightly (>= 34 px total, >= 10 px each side) for 900-dpi or heavily overlapped caps
    for min_tot, min_side in ((42, 12), (34, 10), (26, 8)):
        # Also allow a cap clipped on one side by right/left plot edge or overlapping square
        # if the visible side has >= 22 px and total >= min_tot - 4 (only on outermost cap check)
        cap_mask = (tot_cnt >= min_tot) & (left_cnt >= min_side) & (right_cnt >= min_side)
        if not cap_mask.any():
            continue
        runs = _runs_1d(cap_mask, max_gap=2)
        # Filter out tall solid blobs only if multiple crisp caps exist
        # Actually: if a short error bar has connecting line bridging top & bottom caps into one tall run,
        # the top cap sits at the first ~7 px of the run and bottom cap sits at the last ~7 px!
        top_s, top_e = runs[0]
        bot_s, bot_e = runs[-1]
        if len(runs) >= 2:
            # Measure each cap center from its thin cap band (<= 10 px)
            top_r = rr0 + 0.5 * (top_s + min(top_e, top_s + 9))
            bot_r = rr0 + 0.5 * (max(bot_s, bot_e - 9) + bot_e)
            # Check if a wider cap sits further out with one side partially occluded by another marker
            # (e.g. Fig S2 x=1.01 bottom red cap behind cyan square)
            one_sided = (tot_cnt >= 26) & ((left_cnt >= 17) | (right_cnt >= 17))
            os_runs = _runs_1d(one_sided, max_gap=2)
            if os_runs:
                first_os = os_runs[0]
                last_os = os_runs[-1]
                if first_os[1] - first_os[0] <= 14 and (rr0 + 0.5 * (first_os[0] + first_os[1])) < top_r - 15:
                    top_r = rr0 + 0.5 * (first_os[0] + first_os[1])
                if last_os[1] - last_os[0] <= 14 and (rr0 + 0.5 * (last_os[0] + last_os[1])) > bot_r + 15:
                    bot_r = rr0 + 0.5 * (last_os[0] + last_os[1])
            center = 0.5 * (top_r + bot_r)
            return center, top_r, bot_r, f"bilateral caps (tot>={min_tot})"
        elif top_e - top_s >= 14:
            # Single merged band (short error bar bridged by connecting line or filled disk):
            # take top 6-px cap and bottom 6-px cap of the band!
            top_r = rr0 + 0.5 * (top_s + min(top_e, top_s + 6))
            bot_r = rr0 + 0.5 * (max(top_s, top_e - 6) + top_e)
            center = 0.5 * (top_r + bot_r)
            return center, top_r, bot_r, f"compact bilateral caps (span={top_e - top_s + 1}px)"

    # Fallback: compact marker whose error bar does not protrude past the marker
    stem = mask[rr0:rr1, max(0, col_c - 6) : min(W, col_c + 7)].any(axis=1)
    stem_runs = _runs_1d(stem, max_gap=10)
    if stem_runs:
        best = max(stem_runs, key=lambda t: t[1] - t[0])
        top_r, bot_r = float(rr0 + best[0]), float(rr0 + best[1])
        center = 0.5 * (top_r + bot_r)
        return center, top_r, bot_r, "compact marker/error-bar span"

    raise ValueError(f"No marker/caps detected at col={col_c} in rows [{rr0}, {rr1}]")


def _find_circle_or_square_center_at_col(
    mask: np.ndarray,
    col_c: int,
    r0: int,
    r1: int,
    marker_rad: int = 16,
    filled: bool = False,
) -> float | None:
    """Find the vertical center of a marker (filled disk or open circle/square) at `col_c`."""
    H, W = mask.shape
    c0, c1 = max(0, col_c - marker_rad - 6), min(W, col_c + marker_rad + 7)
    sub = mask[r0:r1, c0:c1]
    if not sub.any():
        return None
    if filled:
        yy, xx = np.ogrid[-marker_rad // 2 : marker_rad // 2 + 1, -marker_rad // 2 : marker_rad // 2 + 1]
        disk = (xx * xx + yy * yy) <= (marker_rad // 2) ** 2
        er = ndimage.binary_erosion(sub, structure=disk)
        rr, cc = np.nonzero(er)
        if len(rr):
            c_rel = col_c - c0
            keep = np.abs(cc - c_rel) <= marker_rad // 2
            if keep.any():
                return float(r0 + rr[keep].mean())
            return float(r0 + rr.mean())
    c_rel = col_c - c0
    left_c = max(0, c_rel - marker_rad)
    right_c = min(sub.shape[1] - 1, c_rel + marker_rad)
    lw = sub[:, max(0, left_c - 3) : min(sub.shape[1], left_c + 4)].any(axis=1)
    rw = sub[:, max(0, right_c - 3) : min(sub.shape[1], right_c + 4)].any(axis=1)
    both = lw & rw
    runs = _runs_1d(both, max_gap=4)
    if runs:
        best = max(runs, key=lambda t: t[1] - t[0])
        return float(r0 + 0.5 * (best[0] + best[1]))
    return None


def calibrate_mathematica_panel(
    img: np.ndarray,
    x_tick_vals: Sequence[float],
    y_tick_vals: Sequence[float],
    x_scale: str = "log10",
    y_scale: str = "linear",
    x_side: str = "bottom",
    y_side: str = "left",
    panel_idx: int | None = None,
) -> tuple[D.Axis, D.Axis, tuple[int, int, int, int], np.ndarray]:
    """Calibrate a Mathematica closed-frame panel with inward ticks.

    For stacked panels like Fig. 5E, pass panel_idx=0 (top) or panel_idx=1 (bottom);
    the x-axis is always calibrated from the bottommost panel's bottom ticks.
    """
    hb, vb, blk = find_frame_spines(img)
    cL, cR = vb[0], vb[-1]
    if panel_idx is None:
        rT, rB = hb[0], hb[-1]
        rB_x = rB
    else:
        rT, rB = hb[panel_idx], hb[panel_idx + 1]
        rB_x = hb[-1]

    xt = find_inward_ticks(blk, hb[0], rB_x, cL, cR, side=x_side, depth=42, min_len=5, max_depth=24)
    xa = match_ticks_regular(xt, x_tick_vals, scale=x_scale)

    yt = find_inward_ticks(blk, rT, rB, cL, cR, side=y_side, depth=42, min_len=5, max_depth=34)
    try:
        ya = match_ticks_regular(yt, y_tick_vals, scale=y_scale)
    except ValueError:
        other = "right" if y_side == "left" else "left"
        yt2 = find_inward_ticks(blk, rT, rB, cL, cR, side=other, depth=42, min_len=5, max_depth=34)
        ya = match_ticks_regular(yt2, y_tick_vals, scale=y_scale)

    box = (cL[1] + 2, rT[1] + 2, cR[0] - 2, rB[0] - 2)
    return xa, ya, box, blk


def read_mathematica_series_at_x(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    box: tuple[int, int, int, int],
    color_mask_arr: np.ndarray,
    x_positions: Sequence[tuple[float, float]],
    legend_rects: Sequence[tuple[int, int, int, int]] = (),
    y_windows: dict[float, tuple[float, float]] | None = None,
    cap_half_w: int = 34,
    filled_rad: int | None = None,
    extra_pad_px: float = 2.0,
    inset_ticks_px: int = 0,
) -> list[dict]:
    """Read a Mathematica marker + symmetric error-bar series at known/detected x-positions.

    `x_positions`: list of (reported_x, plot_axis_x) pairs (equal unless a plot shifts x=0 to tick 10).
    Each marker's true pixel column is refined locally around `xa.to_pixel(plot_axis_x)` from the
    image pixels themselves. Pass `inset_ticks_px=24` for black series so inward x-ticks are excluded.
    """
    c0, r0, c1, r1 = box
    r0_eff = r0 + inset_ticks_px
    r1_eff = r1 - inset_ticks_px
    H, W = img.shape[:2]
    m = np.zeros((H, W), dtype=bool)
    m[r0_eff:r1_eff, c0:c1] = color_mask_arr[r0_eff:r1_eff, c0:c1]
    for lr0, lr1, lc0, lc1 in legend_rects:
        m[max(0, lr0) : min(H, lr1), max(0, lc0) : min(W, lc1)] = False

    rows_out: list[dict] = []
    for rep_x, ax_x in x_positions:
        c_est = int(round(xa.to_pixel(ax_x)))
        y_hint_px = None
        if y_windows and rep_x in y_windows:
            ylo_v, yhi_v = y_windows[rep_x]
            y_hint_px = (max(r0_eff, int(round(ya.to_pixel(yhi_v)))), min(r1_eff, int(round(ya.to_pixel(ylo_v)))))

        rr0 = y_hint_px[0] if y_hint_px else r0_eff
        rr1 = y_hint_px[1] if y_hint_px else r1_eff

        # Refine column center within +/- 4 px using vertical ink density
        win_c0, win_c1 = max(c0, c_est - 4), min(c1, c_est + 5)
        col_prof = m[rr0:rr1, win_c0:win_c1].sum(axis=0)
        if col_prof.max() > 0:
            best_cols = np.nonzero(col_prof >= col_prof.max() * 0.92)[0]
            col_c = int(round(win_c0 + best_cols.mean()))
            center_r, top_r, bot_r, note = _measure_cap_and_center(
                m, col_c, r0_eff, r1_eff, cap_half_w=cap_half_w, y_hint_range=y_hint_px
            )
            half_px = max(2.5, xa.residual_px + ya.residual_px + 2.0)
            # Only refine via filled disk erosion when error bars are shorter than ~45 px
            # (when wide symmetric SEM caps exist, their midpoint is exact and unbiased!)
            if filled_rad is not None and (bot_r - top_r) < 45.0:
                fc = _find_circle_or_square_center_at_col(m, col_c, rr0, rr1, marker_rad=filled_rad, filled=True)
                if fc is not None:
                    center_r = fc
                    note += "; filled disk center"
        else:
            # Marker completely occluded by another curve (e.g. Fig 2C black circle behind red square at x=300):
            col_c = c_est
            left_pts = []
            for dc in range(12, 55, 3):
                cc = col_c - dc
                if cc >= c0:
                    rs = np.nonzero(m[rr0:rr1, cc])[0]
                    if len(rs) and (rs.max() - rs.min()) < 25:
                        left_pts.append((float(cc), float(rr0 + rs.mean())))
            right_pts = []
            for dc in range(12, 55, 3):
                cc = col_c + dc
                if cc < c1:
                    rs = np.nonzero(m[rr0:rr1, cc])[0]
                    if len(rs) and (rs.max() - rs.min()) < 25:
                        right_pts.append((float(cc), float(rr0 + rs.mean())))
            if not left_pts or not right_pts:
                raise ValueError(f"Occluded marker at col={col_c} without flanking connecting line")
            cl, rl = left_pts[0]
            cr, rr_v = right_pts[0]
            t = (col_c - cl) / (cr - cl)
            center_r = rl + t * (rr_v - rl)
            top_r = center_r - 15.0
            bot_r = center_r + 15.0
            half_px = 12.0
            note = "occluded behind overlying marker; interpolated from connecting line"

        x_val = rep_x if rep_x == 0.0 else float(xa.to_value(col_c))
        if rep_x != 0.0 and abs(x_val - rep_x) / abs(rep_x) < 0.025:
            x_val = float(rep_x)

        y_val = ya.to_value(center_r)
        rows_out.append(
            _make_row(
                x=x_val,
                y=y_val,
                col=float(col_c),
                row=float(center_r),
                ya=ya,
                half_px=half_px,
                bar_lo_row=bot_r,
                bar_hi_row=top_r,
                bar_pad_px=extra_pad_px + ya.residual_px,
                note=note,
            )
        )
    return rows_out


# ---------------------------------------------------------------------------
# Specific Matplotlib Panel Helpers (Figs. 3B, 3C, 4B, S3A, S3B, S3C, S4C)
# ---------------------------------------------------------------------------

def calibrate_mpl_panel(
    img: np.ndarray,
    x_tick_vals: Sequence[float] | None,
    y_tick_vals: Sequence[float],
    x_scale: str = "linear",
    y_scale: str = "linear",
    y_side: str = "left",
    x_min_len: int = 8,
    y_min_len: int = 8,
) -> tuple[D.Axis | None, D.Axis, tuple[int, int, int, int], np.ndarray]:
    """Calibrate a Matplotlib closed-frame panel with outward ticks."""
    hb, vb, blk = find_frame_spines(img)
    rT, rB = hb[0], hb[-1]
    cL, cR = vb[0], vb[-1]

    xa = None
    if x_tick_vals is not None:
        xt = find_outward_ticks(blk, rT, rB, cL, cR, side="bottom", depth=65, min_len=x_min_len, max_depth=62)
        xa = match_ticks_regular(xt, x_tick_vals, scale=x_scale)

    yt = find_outward_ticks(blk, rT, rB, cL, cR, side=y_side, depth=65, min_len=y_min_len, max_depth=62)
    ya = match_ticks_regular(yt, y_tick_vals, scale=y_scale)

    box = (cL[1] + 2, rT[1] + 2, cR[0] - 2, rB[0] - 2)
    return xa, ya, box, blk


def read_fig3b_circles(img: np.ndarray, at: Sequence[float]) -> list[dict]:
    """Read the 20 blue circular markers of Fig. 3B (n_g = 1..20 sites)."""
    xa, ya, box, _ = calibrate_mpl_panel(
        img,
        x_tick_vals=[0, 5, 10, 15, 20],
        y_tick_vals=[0.5, 0.4, 0.3, 0.2, 0.1, 0.0],
        x_scale="linear",
        y_scale="linear",
        y_side="left",
    )
    assert xa is not None
    c0, r0, c1, r1 = box
    blue = rgb_dist_mask(img, (31, 119, 180), tol=55.0)
    # Pad blue mask above r0 and below r1? No: markers at y=0 (row 1192) and y=0.49 (row 225)
    # are well inside [r0=186, r1=1232] (circle radius ~59 px; at y=0, bottom hits spine 1234 -> 17 px clip;
    # at n=18..20, top hits spine 184 -> 17 px clip).
    # Measure unclipped circle radius R from n=7..12 using vertical distance transform!
    sub = blue[r0:r1, c0:c1]
    dt = ndimage.distance_transform_edt(sub)
    # Estimate true circle radius from n=7..12
    rads = []
    for n_ref in (7, 8, 9, 10, 11, 12):
        cj = int(round(xa.to_pixel(n_ref))) - c0
        if 0 <= cj < dt.shape[1]:
            rads.append(float(dt[:, max(0, cj - 2) : min(dt.shape[1], cj + 3)].max()))
    R = float(np.median(rads)) if rads else 59.0

    xs = [float(x) for x in at] if at else [float(k) for k in range(1, 21)]
    out: list[dict] = []
    for x in xs:
        col = xa.to_pixel(x)
        cj = int(round(col)) - c0
        cj = max(0, min(sub.shape[1] - 1, cj))
        # Find vertical blue span at this column
        col_slice = sub[:, max(0, cj - 2) : min(sub.shape[1], cj + 3)]
        dt_slice = dt[:, max(0, cj - 2) : min(sub.shape[1], cj + 3)]
        rr_any = np.nonzero(col_slice.any(axis=1))[0]
        if len(rr_any) == 0:
            out.append({"x": x, "y": None, "note": "no blue marker pixels at column"})
            continue

        top_edge = float(rr_any.min())
        bot_edge = float(rr_any.max())
        # Check if clipped by bottom spine (bot_edge >= sub.shape[0] - 3) or top spine (top_edge <= 2)
        if bot_edge >= sub.shape[0] - 3 and (bot_edge - top_edge) < 2.2 * R:
            # Bottom clipped by bottom spine: center is top_edge + R
            rc = r0 + top_edge + R
            note = f"bottom spine clipped; center = top_edge + R ({R:.1f} px)"
            half_px = 4.5
        elif top_edge <= 2 and (bot_edge - top_edge) < 2.2 * R:
            # Top clipped by top spine: center is bot_edge - R
            rc = r0 + bot_edge - R
            note = f"top spine clipped; center = bot_edge - R ({R:.1f} px)"
            half_px = 4.5
        elif col >= c1 - R * 0.8:
            # Rightmost circle (n=20) clipped by right spine: use top/bot edge midpoint!
            rc = r0 + 0.5 * (top_edge + bot_edge)
            if top_edge <= 2:
                rc = r0 + bot_edge - R
            note = "right edge marker; vertical edge fit"
            half_px = 4.5
        else:
            # Interior circle: peak of distance transform gives exact circle center (immune to connecting line!)
            mx = float(dt_slice.max())
            rr_pk = np.nonzero(dt_slice >= mx - 0.5)[0]
            rc = r0 + float(rr_pk.mean())
            note = f"distance-transform circle center (R={mx:.1f} px)"
            half_px = 3.0

        y_val = ya.to_value(rc)
        # Clamp tiny negative values at y=0 (e.g. -0.0005 from 0.5 px) only for note; keep true value
        out.append(_make_row(x=x, y=y_val, col=col, row=rc, ya=ya, half_px=half_px, note=note))
    return out


def read_bar_chart(
    img: np.ndarray,
    y_tick_vals: Sequence[float],
    color: str,
    x_values: Sequence[float],
) -> list[dict]:
    """Read a Matplotlib vertical bar chart (Fig. 4B or Fig. S4C) with centered black SEM error bars."""
    _, ya, box, _ = calibrate_mpl_panel(
        img,
        x_tick_vals=None,
        y_tick_vals=y_tick_vals,
        y_scale="linear",
        y_side="left",
    )
    c0, r0, c1, r1 = box
    r_base = int(round(ya.to_pixel(0.0)))
    r_g, g_g, b_g = _rgb(img)
    if color == "red":
        bar_m = (r_g > 210) & (g_g < 60) & (b_g < 60)
    elif color == "cyan":
        bar_m = (r_g < 60) & (g_g > 210) & (b_g > 210)
    else:
        raise ValueError(f"Unsupported bar color {color!r}")
    blk = (r_g < 60) & (g_g < 60) & (b_g < 60)

    # Identify bars from solid colored regions just above the x-axis baseline
    base_band = bar_m[r_base - 45 : r_base - 10, c0:c1].mean(axis=0) > 0.8
    bar_spans = [(c0 + s, c0 + e) for s, e in _runs_1d(base_band, max_gap=15) if e - s + 1 > 100]
    if len(bar_spans) != len(x_values):
        raise ValueError(f"Expected {len(x_values)} {color} bars, found {len(bar_spans)}: {bar_spans}")

    out: list[dict] = []
    for x_val, (bs, be) in zip(x_values, bar_spans):
        cm = 0.5 * (bs + be)
        ci = int(round(cm))
        # Trace UPWARD continuously from r_base - 10 on the left & right flanks of the bar
        # (completely avoiding both the legend above and the centered black error bar!)
        flank_l = bar_m[: r_base - 5, bs + 20 : ci - 35].mean(axis=1) > 0.8
        flank_r = bar_m[: r_base - 5, ci + 35 : be - 20].mean(axis=1) > 0.8
        flank = flank_l & flank_r
        top_row = r_base - 10
        while top_row > r0 and flank[top_row - 1]:
            top_row -= 1
        y_val = ya.to_value(float(top_row))

        # Locate centered black vertical error bar crossing or touching top_row
        eb_prof = blk[r0 : r_base - 5, ci - 6 : ci + 7].mean(axis=1) > 0.5
        eb_runs = [(r0 + s, r0 + e) for s, e in _runs_1d(eb_prof, max_gap=3)]
        # Keep only the run that overlaps or touches top_row within +/- 15 px
        near = [t for t in eb_runs if t[0] - 15 <= top_row <= t[1] + 15]
        if near:
            eb_top, eb_bot = min(near, key=lambda t: abs(0.5 * (t[0] + t[1]) - top_row))
            # Check if symmetric; if bottom is barely below top_row (or occluded), keep measured ends
            note = f"bar top row={top_row}, black SEM bar rows=[{eb_top}, {eb_bot}]"
        else:
            eb_top, eb_bot = top_row, top_row
            note = f"bar top row={top_row}, SEM error bar <= 1 px (nearly zero)"

        out.append(
            _make_row(
                x=float(x_val),
                y=y_val,
                col=cm,
                row=float(top_row),
                ya=ya,
                half_px=max(2.0, ya.residual_px + 1.5),
                bar_lo_row=float(eb_bot),
                bar_hi_row=float(eb_top),
                bar_pad_px=2.0 + ya.residual_px,
                note=note,
            )
        )
    return out


def read_mpl_curve_at(
    img: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    box: tuple[int, int, int, int],
    color_mask_arr: np.ndarray,
    at: Sequence[float],
    legend_rects: Sequence[tuple[int, int, int, int]] = (),
    half_win_cols: int = 3,
    min_pad_px: float = 2.0,
) -> list[dict]:
    """Read a smooth Matplotlib curve at specified `at` x-values."""
    c0, r0, c1, r1 = box
    H, W = img.shape[:2]
    m = np.zeros((H, W), dtype=bool)
    m[r0:r1, c0:c1] = color_mask_arr[r0:r1, c0:c1]
    for lr0, lr1, lc0, lc1 in legend_rects:
        m[max(0, lr0) : min(H, lr1), max(0, lc0) : min(W, lc1)] = False

    out: list[dict] = []
    for x in at:
        col = xa.to_pixel(float(x))
        ci = int(round(col))
        sub = m[r0:r1, max(c0, ci - half_win_cols) : min(c1, ci + half_win_cols + 1)]
        rr = np.nonzero(sub.any(axis=1))[0]
        if len(rr) == 0:
            # Widen search slightly (+/- 6 px) if a vertical dashed line interrupts the curve at exact column
            sub2 = m[r0:r1, max(c0, ci - 12) : min(c1, ci + 13)]
            rr = np.nonzero(sub2.any(axis=1))[0]
            if len(rr) == 0:
                out.append({"x": float(x), "y": None, "note": "curve not visible at x"})
                continue
            note = "interpolated across narrow gap/line"
            pad = 4.0
        else:
            note = "curve vertical center"
            pad = min_pad_px

        # Pick largest contiguous vertical run (in case stray pixel exists)
        runs = _runs_1d(sub.any(axis=1) if len(sub.nonzero()[0]) else sub2.any(axis=1), max_gap=4)
        best = max(runs, key=lambda t: t[1] - t[0])
        rm = r0 + 0.5 * (best[0] + best[1])
        half_thick = 0.5 * (best[1] - best[0])
        half_px = max(pad, 0.5 * half_thick + ya.residual_px + 1.5)
        y_val = ya.to_value(rm)
        out.append(_make_row(x=float(x), y=y_val, col=col, row=rm, ya=ya, half_px=half_px, note=note))
    return out
