"""Shared image-based calibration and feature-extraction helpers for mcdougal2021 figures.

Every axis coordinate, tick mark, curve height, marker centroid, and error-bar cap
is detected directly from the rendered PDF crop image (no values hard-coded from text).
"""
from __future__ import annotations

import math
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve
from boxgraph.digitize import Axis


def neutral_dark_mask(img: np.ndarray, dark_thr: float = 120.0, max_chroma: float = 35.0) -> np.ndarray:
    """Boolean mask of dark neutral (black/dark-gray) ink pixels."""
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    gray = (r + g + b) / 3.0
    chroma = np.maximum(r, np.maximum(g, b)) - np.minimum(r, np.minimum(g, b))
    return (gray < dark_thr) & (chroma < max_chroma)


def find_box(img: np.ndarray, dark_thr: float = 120.0, neutral_only: bool = True,
             h_min_frac: float = 0.45, v_min_frac: float = 0.75) -> tuple[float, float, float, float]:
    """Locate the rectangular plot frame (c_left, r_top, c_right, r_bot) in continuous coordinates."""
    if neutral_only:
        d = neutral_dark_mask(img, dark_thr=dark_thr)
    else:
        d = img.astype(float).mean(axis=2) < dark_thr
    H, W = d.shape

    row_frac = d.mean(axis=1)
    h_idx = np.nonzero(row_frac > h_min_frac)[0]
    if h_idx.size == 0:
        raise RuntimeError("Could not find horizontal frame borders")
    splits = np.nonzero(np.diff(h_idx) > 6)[0]
    starts = np.r_[0, splits + 1]
    ends = np.r_[splits, h_idx.size - 1]
    hs = [float(np.mean(h_idx[s:e + 1]) + 0.5) for s, e in zip(starts, ends)]
    r_top, r_bot = hs[0], hs[-1]

    sub = d[int(round(r_top)):int(round(r_bot)) + 1, :]
    col_frac = sub.mean(axis=0)
    v_idx = np.nonzero(col_frac > v_min_frac)[0]
    if v_idx.size == 0:
        raise RuntimeError("Could not find vertical frame borders")
    v_splits = np.nonzero(np.diff(v_idx) > 6)[0]
    v_starts = np.r_[0, v_splits + 1]
    v_ends = np.r_[v_splits, v_idx.size - 1]
    vs = [float(np.mean(v_idx[s:e + 1]) + 0.5) for s, e in zip(v_starts, v_ends)]
    return vs[0], r_top, vs[-1], r_bot


def find_ticks_mask(mask: np.ndarray, axis_pos: float, lo: float, hi: float,
                    direction: str, side: int, offset_lo: int = 10, offset_hi: int = 22,
                    min_fill: float = 0.6, max_w: int = 24) -> list[float]:
    """Detect tick centers along an axis line inside a binary ink mask.
    direction='x' scans columns along a horizontal axis at row=axis_pos.
    direction='y' scans rows along a vertical axis at col=axis_pos.
    side=+1 looks at increasing row/col index; side=-1 at decreasing index.
    """
    H, W = mask.shape
    lo_i = max(0, int(round(lo)))
    hi_i = min((W if direction == "x" else H) - 1, int(round(hi)))
    if direction == "x":
        r0 = max(0, min(H - 1, int(round(axis_pos + side * offset_lo))))
        r1 = max(0, min(H - 1, int(round(axis_pos + side * offset_hi))))
        band = mask[min(r0, r1):max(r0, r1) + 1, lo_i:hi_i + 1]
        proj = band.mean(axis=0)
    else:
        c0 = max(0, min(W - 1, int(round(axis_pos + side * offset_lo))))
        c1 = max(0, min(W - 1, int(round(axis_pos + side * offset_hi))))
        band = mask[lo_i:hi_i + 1, min(c0, c1):max(c0, c1) + 1]
        proj = band.mean(axis=1)
    idx = np.nonzero(proj >= min_fill)[0]
    if idx.size == 0:
        return []
    splits = np.nonzero(np.diff(idx) > 4)[0]
    starts = np.r_[0, splits + 1]
    ends = np.r_[splits, idx.size - 1]
    ticks = []
    for s, e in zip(starts, ends):
        w = int(idx[e] - idx[s] + 1)
        if w <= max_w:
            ticks.append(float(np.mean(idx[s:e + 1]) + lo_i + 0.5))
    return ticks


def interval_from_px(axis: Axis, p_center: float, half_px: float) -> tuple[float, float, float]:
    """Convert a pixel position and symmetric pixel half-width to (value, lo, hi) with lo <= value <= hi."""
    val = float(axis.to_value(p_center))
    a = float(axis.to_value(p_center - half_px))
    b = float(axis.to_value(p_center + half_px))
    lo, hi = min(a, b, val), max(a, b, val)
    return val, lo, hi


def trace_black_curve_at(img: np.ndarray, x_axis: Axis, y_axis: Axis,
                         box: tuple[float, float, float, float],
                         xs: list[float],
                         exclude_red: bool = True,
                         y_min_row: int | None = None,
                         y_max_row: int | None = None,
                         tick_pad: int = 28,
                         dark_thr: float = 110.0) -> list[dict]:
    """Trace a single-valued black curve at given x values (in x_axis units).
    Returns BoxGraph mode:code reading dicts with honest 95% intervals.
    """
    c0, r0, c1, r1 = box
    r_lo = int(math.ceil(r0 + 8)) if y_min_row is None else int(y_min_row)
    r_hi = int(math.floor(r1 - 8)) if y_max_row is None else int(y_max_row)
    r_ch = img[:, :, 0].astype(float)
    g_ch = img[:, :, 1].astype(float)
    b_ch = img[:, :, 2].astype(float)
    gray = (r_ch + g_ch + b_ch) / 3.0
    if exclude_red:
        mask = (gray < dark_thr) & ((r_ch - b_ch) < 35.0) & ((r_ch - g_ch) < 35.0)
    else:
        mask = gray < dark_thr

    cal_err_y = max(y_axis.residual_px, 0.5)
    cal_err_x = max(x_axis.residual_px, 0.5)

    # First pass: extract raw center rows and half-thicknesses
    raw = []
    for xv in xs:
        cf = float(x_axis.to_pixel(xv))
        col = int(round(cf))
        if col <= int(math.ceil(c0)) + 2 or col >= int(math.floor(c1)) - 2:
            raw.append((xv, cf, None, None))
            continue
        col_rows = np.nonzero(mask[r_lo:r_hi, col])[0] + r_lo
        if col_rows.size == 0:
            # try +/- 1 col in case of a single-pixel anti-alias gap
            for dc in (-1, 1, -2, 2):
                cc = col + dc
                if int(c0) + 2 < cc < int(c1) - 2:
                    col_rows = np.nonzero(mask[r_lo:r_hi, cc])[0] + r_lo
                    if col_rows.size > 0:
                        break
        if col_rows.size == 0:
            raw.append((xv, cf, None, None))
            continue
        splits = np.nonzero(np.diff(col_rows) > 3)[0]
        starts = np.r_[0, splits + 1]
        ends = np.r_[splits, col_rows.size - 1]
        runs = [(int(col_rows[s]), int(col_rows[e])) for s, e in zip(starts, ends)]
        interior = [rn for rn in runs if rn[0] > r0 + tick_pad and rn[1] < r1 - tick_pad]
        if interior:
            runs = interior
        # pick run closest to previous valid point if multiple, else longest
        if len(runs) == 1:
            ra, rb = runs[0]
        else:
            prev_r = next((p[2] for p in reversed(raw) if p[2] is not None), None)
            if prev_r is not None:
                ra, rb = min(runs, key=lambda rn: abs(0.5 * (rn[0] + rn[1]) - prev_r))
            else:
                ra, rb = max(runs, key=lambda rn: rn[1] - rn[0])
        rmid = 0.5 * (ra + rb) + 0.5
        half_stroke = 0.5 * (rb - ra + 1)
        raw.append((xv, cf, rmid, half_stroke))

    # Second pass: estimate local slope and compute honest 95% intervals
    out = []
    for i, (xv, cf, rmid, half_stroke) in enumerate(raw):
        if rmid is None:
            out.append({"x": float(xv), "y": None, "note": "curve not found at this x"})
            continue
        # Local slope dr/dc
        prev_pt = next(((raw[j][1], raw[j][2]) for j in range(i - 1, -1, -1) if raw[j][2] is not None), None)
        next_pt = next(((raw[j][1], raw[j][2]) for j in range(i + 1, len(raw)) if raw[j][2] is not None), None)
        if prev_pt is not None and next_pt is not None and abs(next_pt[0] - prev_pt[0]) > 1e-6:
            slope = abs((next_pt[1] - prev_pt[1]) / (next_pt[0] - prev_pt[0]))
        elif prev_pt is not None and abs(cf - prev_pt[0]) > 1e-6:
            slope = abs((rmid - prev_pt[1]) / (cf - prev_pt[0]))
        elif next_pt is not None and abs(next_pt[0] - cf) > 1e-6:
            slope = abs((next_pt[1] - rmid) / (next_pt[0] - cf))
        else:
            slope = 0.0
        # Combine half-stroke (clamped conservatively), pixel quantization, calibration misfit, and x-slope propagation
        dx_err = 0.5 + cal_err_x
        half_px = max(half_stroke, 1.5) + 0.5 + cal_err_y + min(slope * dx_err, 12.0)
        yv, lo, hi = interval_from_px(y_axis, rmid, half_px)
        out.append({
            "x": float(xv),
            "y": yv,
            "lo": lo,
            "hi": hi,
            "px": [cf, rmid],
            "px_lo": float(rmid - half_px),
            "px_hi": float(rmid + half_px),
            "note": f"stroke {2 * half_stroke:.1f} px; cal res x={x_axis.residual_px:.2f} y={y_axis.residual_px:.2f} px",
        })
    return out


def make_cross_template(size: int = 85, stroke: float = 14.0,
                        neg_radius: float = 34.0, neg_weight: float = -0.6) -> tuple[np.ndarray, np.ndarray]:
    """Synthetic symmetric 'x' cross template for sub-pixel matched filtering."""
    yy, xx = np.mgrid[:size, :size] - (size - 1) / 2.0
    d1 = np.abs(xx - yy) / math.sqrt(2.0)
    d2 = np.abs(xx + yy) / math.sqrt(2.0)
    on = (d1 <= stroke / 2.0) | (d2 <= stroke / 2.0)
    r2 = xx ** 2 + yy ** 2
    off = (~on) & (r2 <= neg_radius ** 2) & (d1 > stroke / 2.0 + 6.0) & (d2 > stroke / 2.0 + 6.0)
    tpl = np.zeros((size, size), dtype=float)
    tpl[on] = 1.0
    tpl[off] = neg_weight
    return tpl, on.astype(float)


def detect_crosses(mask: np.ndarray, box: tuple[float, float, float, float],
                   size: int = 85, stroke: float = 14.0,
                   thr: float = 0.72, min_dist: int = 20,
                   margin: tuple[int, int, int, int] = (25, 25, 25, 25),
                   legend_box: tuple[int, int, int, int] | None = None,
                   neg_radius: float = 34.0, neg_weight: float = -0.6) -> list[tuple[float, float, float]]:
    """Detect 'x' markers in a binary mask via template correlation + 5x5 peak refinement.
    Returns list of (row, col, score) sorted by (col, row).
    margin = (left, top, right, bottom) inset from box borders.
    """
    c0, r0, c1, r1 = box
    ml, mt, mr, mb = margin
    m = np.zeros_like(mask, dtype=float)
    rs = max(0, int(round(r0)) + mt)
    re = min(mask.shape[0], int(round(r1)) - mb)
    cs = max(0, int(round(c0)) + ml)
    ce = min(mask.shape[1], int(round(c1)) - mr)
    m[rs:re, cs:ce] = mask[rs:re, cs:ce].astype(float)
    if legend_box is not None:
        lc0, lr0, lc1, lr1 = legend_box
        m[lr0:lr1, lc0:lc1] = 0.0

    tpl, on = make_cross_template(size=size, stroke=stroke, neg_radius=neg_radius, neg_weight=neg_weight)
    ideal = float(on.sum())
    sc = fftconvolve(m, tpl[::-1, ::-1], mode="same") / ideal
    loc = ndimage.maximum_filter(sc, size=min_dist)
    pk = (sc == loc) & (sc >= thr)
    rr, cc = np.nonzero(pk)
    order = np.argsort(-sc[rr, cc])
    pts: list[tuple[float, float, float]] = []
    for idx in order:
        ri, ci = int(rr[idx]), int(cc[idx])
        if all((ri - pr) ** 2 + (ci - pc) ** 2 >= min_dist ** 2 for pr, pc, _ in pts):
            # 5x5 score-weighted sub-pixel centroid refinement
            p5 = sc[max(0, ri - 2):ri + 3, max(0, ci - 2):ci + 3]
            w = np.maximum(p5 - p5.min(), 0.0)
            r_ref, c_ref = float(ri), float(ci)
            if w.shape == (5, 5) and w.sum() > 0:
                offs = np.arange(-2, 3, dtype=float)
                r_ref = ri + float((w * offs[:, None]).sum() / w.sum())
                c_ref = ci + float((w * offs[None, :]).sum() / w.sum())
            pts.append((r_ref + 0.5, c_ref + 0.5, float(sc[ri, ci])))
    return sorted(pts, key=lambda p: (p[1], p[0]))
