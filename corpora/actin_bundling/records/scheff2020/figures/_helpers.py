"""Shared fast image-calibration and digitization helpers for scheff2020 figures.

All functions operate strictly on the RGB numpy array `img` passed by BoxGraph
(never re-rendering PDFs or loading external files) and use separable 1D
morphology so every reading executes in < 0.3 s.
"""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage


def clusters_1d(idxs: np.ndarray, max_gap: int = 3) -> list[tuple[int, int, float]]:
    """Group sorted 1D integer indices into contiguous clusters (start, end, center)."""
    if len(idxs) == 0:
        return []
    out = []
    s = p = int(idxs[0])
    for x in idxs[1:]:
        xi = int(x)
        if xi <= p + max_gap:
            p = xi
        else:
            out.append((s, p, 0.5 * (s + p) + 0.5))
            s = p = xi
    out.append((s, p, 0.5 * (s + p) + 0.5))
    return out


def open_rect(mask: np.ndarray, h: int, w: int) -> np.ndarray:
    """Fast separable 2D rectangular morphological opening (h x w)."""
    m = mask.astype(np.uint8)
    if h > 1:
        m = ndimage.minimum_filter1d(m, size=h, axis=0, mode="constant", cval=0)
    if w > 1:
        m = ndimage.minimum_filter1d(m, size=w, axis=1, mode="constant", cval=0)
    if h > 1:
        m = ndimage.maximum_filter1d(m, size=h, axis=0, mode="constant", cval=0)
    if w > 1:
        m = ndimage.maximum_filter1d(m, size=w, axis=1, mode="constant", cval=0)
    return m.astype(bool)


def close_rect(mask: np.ndarray, h: int, w: int) -> np.ndarray:
    """Fast separable 2D rectangular morphological closing (h x w)."""
    m = mask.astype(np.uint8)
    if h > 1:
        m = ndimage.maximum_filter1d(m, size=h, axis=0, mode="constant", cval=0)
    if w > 1:
        m = ndimage.maximum_filter1d(m, size=w, axis=1, mode="constant", cval=0)
    if h > 1:
        m = ndimage.minimum_filter1d(m, size=h, axis=0, mode="constant", cval=0)
    if w > 1:
        m = ndimage.minimum_filter1d(m, size=w, axis=1, mode="constant", cval=0)
    return m.astype(bool)


def find_box(img: np.ndarray, min_frac: float = 0.55, dark_thresh: int = 160) -> dict:
    """Find outer plot frame spines (supports black and colored spines)."""
    ink = img.min(axis=2) < dark_thresh
    rc = clusters_1d(np.where(ink.mean(axis=1) > min_frac)[0])
    cc = clusters_1d(np.where(ink.mean(axis=0) > min_frac)[0])
    if len(rc) < 2 or len(cc) < 2:
        raise ValueError(f"Could not detect 4 plot spines (rows={len(rc)}, cols={len(cc)})")
    return {
        "c0_s": cc[0][0], "c0_e": cc[0][1], "c0": cc[0][2],
        "c1_s": cc[-1][0], "c1_e": cc[-1][1], "c1": cc[-1][2],
        "r0_s": rc[0][0], "r0_e": rc[0][1], "r0": rc[0][2],
        "r1_s": rc[-1][0], "r1_e": rc[-1][1], "r1": rc[-1][2],
    }


def find_ticks(
    img: np.ndarray,
    box: dict,
    side: str,
    outward: bool = False,
    min_len: int = 10,
    max_len: int = 90,
    max_width: int = 24,
    corner_pad: int = 12,
    dark_thresh: int = 160,
) -> list[tuple[float, int]]:
    """Detect tick positions and lengths along a spine ('bottom', 'top', 'left', 'right')."""
    ink = img.min(axis=2) < dark_thresh
    H, W = ink.shape
    if side in ("bottom", "top"):
        c_lo = box["c0_e"] + corner_pad
        c_hi = box["c1_s"] - corner_pad
        if side == "bottom":
            r_start = box["r1_e"] + 1 if outward else box["r1_s"] - 1
            step = 1 if outward else -1
        else:
            r_start = box["r0_s"] - 1 if outward else box["r0_e"] + 1
            step = -1 if outward else 1
        # Vectorized contiguous run length along rows
        max_k = min(max_len, (H - r_start) if step > 0 else (r_start + 1))
        rows = r_start + step * np.arange(max_k)
        strip = ink[rows, c_lo:c_hi]  # (max_k, N)
        # First False along axis 0 gives run length
        has_false = ~strip
        first_false = np.argmax(has_false, axis=0)
        all_true = strip.all(axis=0)
        lens = np.where(all_true, max_k, first_false)
        hits = np.where((lens >= min_len) & (lens < max_len - 4))[0]
        res = []
        for s, e, mid in clusters_1d(hits, max_gap=2):
            if e - s + 1 <= max_width:
                res.append((c_lo + mid, int(lens[s:e + 1].max())))
        return res
    else:
        r_lo = box["r0_e"] + corner_pad
        r_hi = box["r1_s"] - corner_pad
        if side == "left":
            c_start = box["c0_s"] - 1 if outward else box["c0_e"] + 1
            step = -1 if outward else 1
        else:
            c_start = box["c1_e"] + 1 if outward else box["c1_s"] - 1
            step = 1 if outward else -1
        max_k = min(max_len, (W - c_start) if step > 0 else (c_start + 1))
        cols = c_start + step * np.arange(max_k)
        strip = ink[r_lo:r_hi, cols]  # (N, max_k)
        has_false = ~strip
        first_false = np.argmax(has_false, axis=1)
        all_true = strip.all(axis=1)
        lens = np.where(all_true, max_k, first_false)
        hits = np.where((lens >= min_len) & (lens < max_len - 4))[0]
        res = []
        for s, e, mid in clusters_1d(hits, max_gap=2):
            if e - s + 1 <= max_width:
                res.append((r_lo + mid, int(lens[s:e + 1].max())))
        return res


class AxisCal:
    """Linear or log10 axis calibration from (pixel, value) pairs."""

    def __init__(self, pairs: list[tuple[float, float]], scale: str = "linear", max_res_px: float = 3.0):
        self.scale = scale
        self.pairs = [(float(p), float(v)) for p, v in pairs]
        ps = np.array([p for p, _ in self.pairs], dtype=float)
        if scale == "log10":
            us = np.array([math.log10(v) for _, v in self.pairs], dtype=float)
        else:
            us = np.array([v for _, v in self.pairs], dtype=float)
        self.b, self.a = np.polyfit(ps, us, 1)  # u = a + b * p
        self.residual_px = float(np.max(np.abs(ps - (us - self.a) / self.b)))
        if self.residual_px > max_res_px:
            raise ValueError(f"Axis calibration residual {self.residual_px:.2f} px exceeds {max_res_px} px")

    def to_value(self, p: float) -> float:
        u = self.a + self.b * float(p)
        return float(10.0 ** u) if self.scale == "log10" else float(u)

    def to_pixel(self, v: float) -> float:
        u = math.log10(float(v)) if self.scale == "log10" else float(v)
        return float((u - self.a) / self.b)


def find_inner_frame_box(img: np.ndarray, box: dict, min_w_frac: float = 0.15, min_h_frac: float = 0.10) -> tuple[int, int, int, int] | None:
    """Detect a black-bordered legend rectangle inside the plot box, returning (c0, r0, c1, r1) or None."""
    r0i, r1i = box["r0_e"] + 15, box["r1_s"] - 15
    c0i, c1i = box["c0_e"] + 15, box["c1_s"] - 15
    # 4x downsampling makes legend border detection < 5 ms while retaining ~4 px accuracy
    step = 4
    sub = img[r0i:r1i:step, c0i:c1i:step]
    dark = (sub.max(axis=2) < 110)
    H, W = dark.shape
    min_w = max(8, int(min_w_frac * W))
    min_h = max(8, int(min_h_frac * H))
    h_opened = open_rect(dark, 1, min_w)
    v_opened = open_rect(dark, min_h, 1)
    frame = h_opened | v_opened
    lab, n = ndimage.label(frame)
    best = None
    best_area = 0
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        h = rr.max() - rr.min() + 1
        w = cc.max() - cc.min() + 1
        if w >= min_w and h >= min_h and w * h > best_area:
            best_area = w * h
            best = (
                c0i + int(cc.min()) * step - 8,
                r0i + int(rr.min()) * step - 8,
                c0i + int(cc.max() + 1) * step + 8,
                r0i + int(rr.max() + 1) * step + 8,
            )
    return best


def labeled_components(mask: np.ndarray):
    """Yield (label_id, r_slice, c_slice, local_bool_mask) using fast ndimage.find_objects."""
    lab, n = ndimage.label(mask)
    objs = ndimage.find_objects(lab)
    for i, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        r_sl, c_sl = sl
        loc = (lab[r_sl, c_sl] == i)
        yield i, r_sl, c_sl, loc


def sample_curve_at(
    mask: np.ndarray,
    xa: AxisCal,
    ya: AxisCal,
    at: list[float],
    box: dict,
    fit_half_win: int = 55,
    max_gap_px: int = 120,
    min_pts: int = 5,
) -> list[dict]:
    """Sample a single-valued solid or dashed curve at requested x positions.

    Uses per-column median rows of `mask` inside `box`. When ink is present
    within ±3 px of `c_at`, reads directly from the local band; when `c_at`
    falls in a dash gap or behind an occluding marker, fits a local weighted
    linear (or quadratic if wide) regression to neighboring ink columns and
    widens the 95% reading interval honestly.
    """
    r_lo = box["r0_e"] + 4
    r_hi = box["r1_s"] - 4
    c_lo = box["c0_e"] + 4
    c_hi = box["c1_s"] - 4

    cols_list = []
    rows_med = []
    half_thk = []
    for c in range(c_lo, c_hi):
        rs = np.nonzero(mask[r_lo:r_hi, c])[0]
        if len(rs) == 0:
            continue
        rs = rs + r_lo
        # Reject outlier stray specks far from median if multiple runs exist
        runs = clusters_1d(rs, max_gap=2)
        if len(runs) == 1:
            s, e, mid = runs[0]
            if e - s + 1 <= 45:
                cols_list.append(c + 0.5)
                rows_med.append(mid)
                half_thk.append(0.5 * (e - s + 1))
        else:
            # Keep runs <= 45 px; pick the longest run (most likely curve ink)
            good = [r for r in runs if (r[1] - r[0] + 1) <= 45]
            if good:
                best = max(good, key=lambda r: r[1] - r[0] + 1)
                cols_list.append(c + 0.5)
                rows_med.append(best[2])
                half_thk.append(0.5 * (best[1] - best[0] + 1))

    cols_arr = np.array(cols_list, dtype=float)
    rows_arr = np.array(rows_med, dtype=float)
    thk_arr = np.array(half_thk, dtype=float)
    med_thk = float(np.median(thk_arr)) if len(thk_arr) else 3.0

    # Light median filter pass to reject single-column stray specks
    if len(rows_arr) >= 7:
        # Compare each point to median of nearby points within +-25 px
        keep = np.ones(len(rows_arr), dtype=bool)
        for i in range(len(rows_arr)):
            near = np.abs(cols_arr - cols_arr[i]) <= 30.0
            if near.sum() >= 5:
                if abs(rows_arr[i] - np.median(rows_arr[near])) > max(25.0, 4.0 * med_thk):
                    keep[i] = False
        cols_arr, rows_arr, thk_arr = cols_arr[keep], rows_arr[keep], thk_arr[keep]

    readings = []
    for x in at:
        c_at = xa.to_pixel(x)
        if len(cols_arr) == 0:
            readings.append({"x": float(x), "y": None, "note": "no curve pixels found"})
            continue
        dist = np.abs(cols_arr - c_at)
        min_d = float(dist.min())
        if min_d > max_gap_px:
            readings.append({"x": float(x), "y": None, "note": f"no curve ink within {max_gap_px} px"})
            continue

        # Direct hit within +-3 px
        direct = dist <= 3.0
        if direct.sum() >= 2:
            r_est = float(np.median(rows_arr[direct]))
            loc_thk = float(np.median(thk_arr[direct]))
            # Estimate local slope for x-calibration propagation
            win = dist <= fit_half_win
            slope = 0.0
            if win.sum() >= min_pts and np.ptp(cols_arr[win]) > 6.0:
                slope = float(np.polyfit(cols_arr[win] - c_at, rows_arr[win], 1)[0])
            dh_px = loc_thk + 0.5 + ya.residual_px + abs(slope) * (0.5 + xa.residual_px)
            note = f"direct ink ({int(direct.sum())} cols, half_w={loc_thk:.1f} px)"
        else:
            # Interpolate across dash gap / marker occlusion using weighted local linear/quad fit
            win = dist <= fit_half_win
            if win.sum() < min_pts or np.ptp(cols_arr[win]) < 8.0:
                win = dist <= max_gap_px
            if win.sum() < 2:
                readings.append({"x": float(x), "y": None, "note": "too few nearby ink columns"})
                continue
            dx = cols_arr[win] - c_at
            ys = rows_arr[win]
            # Gaussian weights favoring closest columns on both sides
            sig = max(25.0, min_d * 1.2)
            w = np.exp(-0.5 * (dx / sig) ** 2)
            # Use degree 2 if we have points on both sides and >= 12 cols, else degree 1
            has_left = (dx < 0).sum() >= 4
            has_right = (dx > 0).sum() >= 4
            deg = 2 if (has_left and has_right and win.sum() >= 12) else 1
            coeffs = np.polyfit(dx, ys, deg, w=np.sqrt(w))
            r_est = float(coeffs[-1])
            slope = float(coeffs[-2]) if deg >= 1 else 0.0
            fit_res = float(np.sqrt(np.average((np.polyval(coeffs, dx) - ys) ** 2, weights=w)))
            loc_thk = float(np.median(thk_arr[win]))
            gap_pen = 0.08 * min_d + 2.0 * fit_res
            dh_px = loc_thk + 0.5 + ya.residual_px + abs(slope) * (0.5 + xa.residual_px) + gap_pen
            note = f"interpolated across {min_d:.0f} px gap (deg={deg}, rms={fit_res:.1f} px)"

        y_val = ya.to_value(r_est)
        y_a = ya.to_value(r_est - dh_px)
        y_b = ya.to_value(r_est + dh_px)
        lo, hi = min(y_a, y_b), max(y_a, y_b)
        readings.append({
            "x": float(x),
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [c_at, r_est],
            "px_lo": r_est - dh_px,
            "px_hi": r_est + dh_px,
            "note": note,
        })
    return readings


def trace_vertical_ebar(
    stem_mask: np.ndarray,
    c_center: float,
    r_top_marker: int,
    r_bot_marker: int,
    r_min_bound: int,
    r_max_bound: int,
    col_half_win: int = 5,
    max_gap: int = 10,
) -> tuple[float | None, float | None]:
    """Trace vertical error bar stem upward from r_top_marker and downward from r_bot_marker.

    Returns (r_upper_end, r_lower_end) in continuous pixel rows, or None if no stem exists
    on that side (or hidden inside marker).
    """
    H, W = stem_mask.shape
    c0 = max(0, int(round(c_center)) - col_half_win)
    c1 = min(W, int(round(c_center)) + col_half_win + 1)

    # For each row, check if there is stem ink near c_center (narrow vertical stroke or horizontal cap)
    # A stem pixel exists at row r if any pixel in [c0, c1) is True
    col_profile = stem_mask[:, c0:c1].any(axis=1)

    # Trace UPWARD from r_top_marker - 1 down to r_min_bound
    r_up = None
    gap = 0
    for r in range(r_top_marker - 1, max(r_min_bound, 0) - 1, -1):
        if col_profile[r]:
            r_up = r + 0.5
            gap = 0
        else:
            gap += 1
            if gap > max_gap:
                break

    # Trace DOWNWARD from r_bot_marker + 1 up to r_max_bound
    r_dn = None
    gap = 0
    for r in range(r_bot_marker + 1, min(r_max_bound, H - 1) + 1):
        if col_profile[r]:
            r_dn = r + 0.5
            gap = 0
        else:
            gap += 1
            if gap > max_gap:
                break

    return r_up, r_dn


def make_marker_reading(
    xa: AxisCal,
    ya: AxisCal,
    c: float,
    r: float,
    r_up: float | None,
    r_dn: float | None,
    marker_half_h: float,
    note: str = "",
    force_x: float | None = None,
    include_ebar: bool = True,
) -> dict:
    """Build a BoxGraph reading dict for a marker, optionally with vertical error bar ends."""
    x_val = float(force_x) if force_x is not None else xa.to_value(c)
    y_val = ya.to_value(r)
    h_px = max(1.5, 0.25 * marker_half_h) + ya.residual_px
    ya1, ya2 = ya.to_value(r - h_px), ya.to_value(r + h_px)
    lo, hi = min(ya1, ya2), max(ya1, ya2)
    out = {
        "x": x_val,
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [float(c), float(r)],
        "px_lo": float(r - h_px),
        "px_hi": float(r + h_px),
        "note": note,
    }
    if include_ebar and (r_up is not None or r_dn is not None):
        # If one side is missing (e.g. clipped at spine or zero upper error bar), handle honestly
        e_up = r_up if r_up is not None else r
        e_dn = r_dn if r_dn is not None else r
        bar_pad = 2.0 + ya.residual_px
        b1, b2 = ya.to_value(e_up), ya.to_value(e_dn)
        b1_out, b2_out = ya.to_value(e_up - bar_pad), ya.to_value(e_dn + bar_pad)
        out["bar_lo"] = min(b1, b2, y_val)
        out["bar_hi"] = max(b1, b2, y_val)
        out["bar_lo_outer"] = min(b1_out, b2_out, lo)
        out["bar_hi_outer"] = max(b1_out, b2_out, hi)
        out["px_lo"] = float(e_up)
        out["px_hi"] = float(e_dn)
    return out
