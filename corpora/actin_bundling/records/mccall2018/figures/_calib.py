"""Shared image-based axis calibration and figure reading helpers for mccall2018.

Every calibration and data reading is derived strictly from the passed RGB `image`
array (no PDF re-rendering inside `read(image, at)` and no hardcoded y values).
"""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage

from boxgraph import digitize as D


# ---------------------------------------------------------------------------
# Generic frame & tick detection utilities
# ---------------------------------------------------------------------------

def _gray_sat(img: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    f = img.astype(np.float32)
    gray = f.mean(axis=2)
    sat = f.max(axis=2) - f.min(axis=2)
    return gray, sat


def _cluster_1d(indices: np.ndarray, gap: int = 8) -> list[tuple[int, int, float]]:
    if len(indices) == 0:
        return []
    out: list[tuple[int, int, float]] = []
    s = int(indices[0])
    p = s
    for x in indices[1:]:
        xi = int(x)
        if xi - p > gap:
            out.append((s, p, 0.5 * (s + p)))
            s = xi
        p = xi
    out.append((s, p, 0.5 * (s + p)))
    return out


def find_frame(
    img: np.ndarray,
    *,
    dark_thr: float = 150.0,
    sat_thr: float = 40.0,
    min_h_frac: float = 0.45,
    min_v_frac: float = 0.45,
    h_region: tuple[float, float] = (0.0, 1.0),
    v_region: tuple[float, float] = (0.0, 1.0),
) -> tuple[float, float, float, float]:
    """Detect rectangular plot frame (cL, rT, cR, rB) using vectorized dark-fraction projections."""
    gray, sat = _gray_sat(img)
    dark = (gray < dark_thr) & (sat < sat_thr)
    H, W = dark.shape

    r0, r1 = int(H * h_region[0]), int(H * h_region[1])
    c0, c1 = int(W * v_region[0]), int(W * v_region[1])

    # Horizontal spines: rows with high dark fraction across [c0:c1]
    row_frac = dark[r0:r1, c0:c1].mean(axis=1)
    h_idx = np.where(row_frac >= min_h_frac)[0] + r0
    h_cls = _cluster_1d(h_idx, gap=8)
    if len(h_cls) < 2:
        raise RuntimeError(f"Expected >=2 horizontal frame spines, found {h_cls}")
    # Keep spines close to the maximum horizontal coverage (excludes short legend borders)
    h_scores = [float(row_frac[a - r0 : b - r0 + 1].max()) for a, b, _ in h_cls]
    max_h = max(h_scores)
    h_good = [c for c, sc in zip(h_cls, h_scores) if sc >= 0.78 * max_h]
    rT = float(h_good[0][2])
    rB = float(h_good[-1][2])

    # Vertical spines restricted to rows [rT..rB]
    ra = max(0, int(round(rT)) - 4)
    rb = min(H, int(round(rB)) + 5)
    col_frac = dark[ra:rb, c0:c1].mean(axis=0)
    v_idx = np.where(col_frac >= max(min_v_frac, 0.55))[0] + c0
    v_cls = _cluster_1d(v_idx, gap=8)
    if len(v_cls) < 2:
        raise RuntimeError(f"Expected >=2 vertical frame spines, found {v_cls}")
    v_scores = [float(col_frac[a - c0 : b - c0 + 1].max()) for a, b, _ in v_cls]
    max_v = max(v_scores)
    v_good = [c for c, sc in zip(v_cls, v_scores) if sc >= 0.78 * max_v]
    cL = float(v_good[0][2])
    cR = float(v_good[-1][2])
    return cL, rT, cR, rB


def find_inward_ticks(
    img: np.ndarray,
    box: tuple[float, float, float, float],
    side: str,
    *,
    dark_thr: float = 150.0,
    sat_thr: float = 40.0,
    depth_start: int = 7,
    depth_end: int = 28,
    min_depth: int = 10,
    max_width: int = 14,
    margin: int = 8,
) -> list[tuple[float, int]]:
    """Detect inward tick marks on 'bottom', 'top', 'left', or 'right' spine.

    Returns list of (center_pixel, tick_depth_px).
    """
    gray, sat = _gray_sat(img)
    dark = (gray < dark_thr) & (sat < sat_thr)
    H, W = dark.shape
    cL, rT, cR, rB = box

    if side == "bottom":
        rs = max(0, int(round(rB)) - depth_end)
        re = max(0, int(round(rB)) - depth_start)
        cs = max(0, int(round(cL)) + margin)
        ce = min(W, int(round(cR)) - margin)
        strip = dark[rs:re, cs:ce]
        offset = cs
        cnt = strip.sum(axis=0)
    elif side == "top":
        rs = min(H, int(round(rT)) + depth_start)
        re = min(H, int(round(rT)) + depth_end)
        cs = max(0, int(round(cL)) + margin)
        ce = min(W, int(round(cR)) - margin)
        strip = dark[rs:re, cs:ce]
        offset = cs
        cnt = strip.sum(axis=0)
    elif side == "left":
        rs = max(0, int(round(rT)) + margin)
        re = min(H, int(round(rB)) - margin)
        cs = min(W, int(round(cL)) + depth_start)
        ce = min(W, int(round(cL)) + depth_end)
        strip = dark[rs:re, cs:ce].T
        offset = rs
        cnt = strip.sum(axis=0)
    elif side == "right":
        rs = max(0, int(round(rT)) + margin)
        re = min(H, int(round(rB)) - margin)
        cs = max(0, int(round(cR)) - depth_end)
        ce = max(0, int(round(cR)) - depth_start)
        strip = dark[rs:re, cs:ce].T
        offset = rs
        cnt = strip.sum(axis=0)
    else:
        raise ValueError(f"Unknown side {side!r}")

    idx = np.where(cnt >= min_depth)[0]
    clusters = _cluster_1d(idx, gap=3)
    out = []
    for a, b, m in clusters:
        w = b - a + 1
        if w <= max_width:
            out.append((offset + m, int(cnt[a : b + 1].max())))
    return out


def match_ticks_near(
    cands: list[tuple[float, int]],
    expected_px: list[float],
    values: list[float],
    *,
    tol_px: float = 15.0,
    min_depth: int = 0,
    scale: str = "linear",
) -> D.Axis:
    """Match detected tick candidates to expected pixel positions and fit Axis."""
    marks: list[tuple[float, float]] = []
    filt = [p for p, d in cands if d >= min_depth]
    for ep, val in zip(expected_px, values):
        near = [p for p in filt if abs(p - ep) <= tol_px]
        if near:
            best = min(near, key=lambda p: abs(p - ep))
            marks.append((best, val))
    if len(marks) < 2:
        raise RuntimeError(f"Too few matched ticks ({len(marks)}): cands={cands[:12]}")
    return D.Axis(scale, marks)


def _interval_from_px(
    ya: D.Axis,
    row: float,
    half_px: float,
) -> tuple[float, float, float, float, float]:
    """Return (y, lo, hi, px_lo, px_hi) from center row and half-width in pixels."""
    y = ya.to_value(row)
    a = ya.to_value(row - half_px)
    b = ya.to_value(row + half_px)
    lo, hi = (min(a, b), max(a, b))
    # Guarantee lo <= y <= hi numerically
    lo = min(lo, y)
    hi = max(hi, y)
    return y, lo, hi, float(row - half_px), float(row + half_px)


def _trace_mask_at(
    mask: np.ndarray,
    xa: D.Axis,
    ya: D.Axis,
    box: tuple[float, float, float, float],
    at: list[float],
    *,
    half_w: int = 3,
    fallback_w: int = 10,
    inset: int = 4,
    extra_px: float = 1.5,
    mode: str = "median",
    note_prefix: str = "",
) -> list[dict]:
    """Read y(x) at requested `at` positions from a boolean curve mask."""
    cL, rT, cR, rB = box
    H, W = mask.shape
    r0 = max(0, int(math.ceil(rT + inset)))
    r1 = min(H, int(math.floor(rB - inset)))
    c0 = max(0, int(math.ceil(cL + inset)))
    c1 = min(W, int(math.floor(cR - inset)))

    out: list[dict] = []
    for xv in at:
        col = xa.to_pixel(xv)
        ci = int(round(col))
        found_rows = np.array([], dtype=float)
        used_w = half_w
        for hw in (half_w, fallback_w):
            a = max(c0, ci - hw)
            b = min(c1, ci + hw + 1)
            if b <= a:
                continue
            ys, xs = np.where(mask[r0:r1, a:b])
            if len(ys) > 0:
                found_rows = ys.astype(float) + r0
                used_w = hw
                break
        if len(found_rows) == 0:
            out.append({"x": float(xv), "y": None, "found": False, "note": "no ink found"})
            continue

        # If multiple vertical clusters exist (e.g., noise speckle), keep largest cluster
        sr = np.sort(found_rows)
        if sr[-1] - sr[0] > 35:
            # cluster gaps > 12 px and pick largest cluster
            gaps = np.where(np.diff(sr) > 12)[0]
            starts = np.r_[0, gaps + 1]
            ends = np.r_[gaps + 1, len(sr)]
            best_i = int(np.argmax(ends - starts))
            sr = sr[starts[best_i] : ends[best_i]]

        if mode == "median":
            row = float(np.median(sr))
        elif mode == "mean":
            row = float(np.mean(sr))
        else:
            row = 0.5 * (float(sr.min()) + float(sr.max()))

        thickness_half = 0.5 * (float(sr.max()) - float(sr.min()))
        cal_err = getattr(ya, "residual_px", 0.5)
        h_px = max(1.5, thickness_half + cal_err + extra_px + (used_w - half_w) * 0.3)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, row, h_px)
        note = f"{note_prefix}cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
        out.append({
            "x": float(xv),
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(col), float(row)],
            "px_lo": float(px_lo),
            "px_hi": float(px_hi),
            "note": note,
        })
    return out


# ---------------------------------------------------------------------------
# Color helpers
# ---------------------------------------------------------------------------

def _rgb_dist_mask(img: np.ndarray, rgb: tuple[float, float, float], tol: float) -> np.ndarray:
    f = img.astype(np.float32)
    t = np.array(rgb, dtype=np.float32)
    d2 = ((f - t) ** 2).sum(axis=2)
    return d2 <= (tol * tol)


# ---------------------------------------------------------------------------
# Figure 1D (page 2, 1200 dpi): Normalized intensity line scans (Actin, BSA, DIC)
# ---------------------------------------------------------------------------

def calibrate_fig1d(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    # X ticks: -6 at cL, +6 at cR, major inward ticks at -4, -2, 0, +2, +4
    b_ticks = find_inward_ticks(img, box, "bottom", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    exp_x = [cL + w * (v + 6.0) / 12.0 for v in (-4.0, -2.0, 0.0, 2.0, 4.0)]
    x_marks = [(cL, -6.0), (cR, 6.0)]
    for ep, val in zip(exp_x, (-4.0, -2.0, 0.0, 2.0, 4.0)):
        near = [p for p, _ in b_ticks if abs(p - ep) <= 15]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y ticks: major inward ticks at 1.0 (~rT+20), 0.5 (~mid), 0.0 (~rB-20), plus 0.1 minor ticks
    l_ticks = find_inward_ticks(img, box, "left", depth_start=6, depth_end=30, min_depth=10, max_width=14)
    # Identify major ticks (depth >= 20)
    maj_y = sorted([p for p, d in l_ticks if d >= 20 and rT + 10 < p < rB - 10])
    # Expect 3 major ticks: 1.0, 0.5, 0.0
    if len(maj_y) >= 3:
        y_marks = [(maj_y[0], 1.0), (maj_y[len(maj_y) // 2], 0.5), (maj_y[-1], 0.0)]
        # Also include clean minor ticks at 0.9, 0.8, ..., 0.1
        ya0 = D.Axis("linear", y_marks)
        for v in (0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1):
            ep = ya0.to_pixel(v)
            near = [p for p, _ in l_ticks if abs(p - ep) <= 10]
            if near:
                y_marks.append((min(near, key=lambda p: abs(p - ep)), v))
        ya = D.Axis("linear", sorted(y_marks))
    else:
        raise RuntimeError(f"Fig 1D major y-ticks not found: {l_ticks}")
    return xa, ya, box


def read_fig1d(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_fig1d(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "actin":
        mask = (g - r > 55) & (g - b > 20) & (g > 100)
    elif series == "bsa":
        mask = (r - b > 90) & (r > 180) & (g > 100) & (g < 220)
    elif series == "dic":
        mask = (gray < 110) & (sat < 35)
    else:
        raise ValueError(series)

    # Mask out legend box (center-top interior: x in [-1.5, 1.4], y in [0.43, 0.79])
    lc0 = int(round(xa.to_pixel(-1.55)))
    lc1 = int(round(xa.to_pixel(1.45)))
    lr0 = int(round(ya.to_pixel(0.80)))
    lr1 = int(round(ya.to_pixel(0.42)))
    mask[lr0:lr1, lc0:lc1] = False

    # Also mask out tick marks near bottom/top/left/right spines for black DIC curve
    if series == "dic":
        mask[int(rB) - 45 :, :] = False
        mask[: int(rT) + 45, :] = False
        mask[:, : int(cL) + 45] = False
        mask[:, int(cR) - 45 :] = False

    return _trace_mask_at(mask, xa, ya, box, at, half_w=4, fallback_w=12, inset=4, extra_px=1.5)


# ---------------------------------------------------------------------------
# Figure 1E (page 2, 1200 dpi): Average RAM bar chart + SEM error bars
# ---------------------------------------------------------------------------

def calibrate_fig1e(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40)
    cL, rT, cR, rB = box
    h = rB - rT
    r_ticks = find_inward_ticks(img, box, "right", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    # Major ticks at 40 (rT), 30, 20, 10, 0 (rB)
    y_marks = [(rT, 40.0), (rB, 0.0)]
    for val in (30.0, 20.0, 10.0):
        ep = rT + h * (40.0 - val) / 40.0
        near = [p for p, _ in r_ticks if abs(p - ep) <= 15]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    ya = D.Axis("linear", sorted(y_marks))
    return ya, box


def read_fig1e(img: np.ndarray, group: str) -> list[dict]:
    ya, box = calibrate_fig1e(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)
    dark = (gray < 120) & (sat < 40)

    if group == "actin":
        cmask = (g - r > 55) & (g - b > 20)
        x_indices = [1.0, 2.0]
        labels = ["TMR-actin + BSA (n=171)", "TMR-actin alone (n=117)"]
    elif group == "bsa":
        cmask = (r - b > 90) & (r > 180) & (g > 100)
        x_indices = [3.0, 4.0]
        labels = ["647-BSA + actin (n=214)", "647-BSA alone (n=287)"]
    else:
        raise ValueError(group)

    cols = np.where(cmask[int(rT) + 10 : int(rB) - 10, int(cL) + 10 : int(cR) - 10].sum(axis=0) > 40)[0]
    cols = cols + int(cL) + 10
    clusters = _cluster_1d(cols, gap=20)
    if len(clusters) != 2:
        raise RuntimeError(f"Fig 1E {group}: expected 2 bars, found {len(clusters)}")

    out: list[dict] = []
    for (x0, x1, xm), xi, lbl in zip(clusters, x_indices, labels):
        bw = x1 - x0
        # Find bar top border row from side columns (away from center error bar stem)
        sub = cmask[int(rT) + 5 : int(rB) - 5, x0 + 4 : x1 - 4]
        mid = sub.shape[1] // 2
        side = np.r_[np.arange(0, max(1, mid - 20)), np.arange(min(sub.shape[1], mid + 20), sub.shape[1])]
        tops = [int(np.where(sub[:, k])[0][0]) + int(rT) + 5 for k in side if sub[:, k].any()]
        color_top = float(np.median(tops))

        # Locate exact center of the horizontal black bar-top stroke right above color_top
        s0 = max(int(rT) + 3, int(round(color_top)) - 12)
        s1 = min(int(rB) - 3, int(round(color_top)) + 4)
        full_w = dark[s0:s1, x0 + 2 : x1 - 2].sum(axis=1)
        border_rows = np.where(full_w > 0.70 * bw)[0] + s0
        bar_top = float(np.mean(border_rows)) if len(border_rows) else color_top

        # Search for SEM error bar upper cap above bar_top (within 35 px) and lower cap below bar_top
        # Exclude full-width bar border and exclude interior n=... text (> 24 px below bar top for orange)
        search_up = 35 if group == "actin" else 14
        search_dn = 35 if group == "actin" else 12
        e0 = max(int(rT) + 3, int(round(bar_top)) - search_up)
        e1 = min(int(rB) - 3, int(round(bar_top)) + search_dn)

        # Look for dark pixels in cap band [xm-25 .. xm+25] outside full-width bar border
        cap_counts = dark[e0:e1, int(round(xm)) - 25 : int(round(xm)) + 26].sum(axis=1)
        full_counts = dark[e0:e1, x0 + 2 : x1 - 2].sum(axis=1)
        # Upper cap: rows above bar_top (before full_counts hits > 0.65*bw) with cap_counts >= 14
        upper_cands = []
        lower_cands = []
        for idx_r in range(e1 - e0):
            rr = e0 + idx_r
            if full_counts[idx_r] > 0.65 * bw:
                continue
            if cap_counts[idx_r] >= 14:
                if rr < bar_top - 1.5:
                    upper_cands.append(rr)
                elif rr > bar_top + 1.5:
                    lower_cands.append(rr)

        if upper_cands:
            # Take topmost cluster center for upper cap
            uc = _cluster_1d(np.array(upper_cands), gap=2)[0][2]
        else:
            # Cap touches top of bar border (tiny SEM): take topmost dark pixel on center stem
            stem = np.where(dark[e0 : int(round(bar_top)) + 1, int(round(xm)) - 10 : int(round(xm)) + 11].sum(axis=1) >= 8)[0] + e0
            uc = float(stem[0]) if len(stem) else bar_top - 2.0

        if lower_cands:
            lc = _cluster_1d(np.array(lower_cands), gap=2)[-1][2]
        else:
            # Lower cap hidden/merged with bar top stroke (orange bars): mirror upper cap
            lc = 2.0 * bar_top - uc

        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, bar_top, max(1.5, ya.residual_px + 1.5))
        bar_hi = ya.to_value(uc)
        bar_lo = ya.to_value(lc)
        end_err_px = max(1.5, ya.residual_px + 1.5)
        bar_hi_outer = ya.to_value(uc - end_err_px)
        bar_lo_outer = ya.to_value(lc + end_err_px)
        bar_lo, bar_hi = min(bar_lo, y), max(bar_hi, y)
        bar_lo_outer = min(bar_lo_outer, bar_lo)
        bar_hi_outer = max(bar_hi_outer, bar_hi)
        sem_note = f"{lbl}; bar_top={bar_top:.1f}, caps=({uc:.1f},{lc:.1f}), cal_res={ya.residual_px:.2f}px"
        if not lower_cands:
            sem_note += " (lower cap mirrored across bar top stroke)"
        out.append({
            "x": float(xi),
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "bar_lo": float(bar_lo),
            "bar_hi": float(bar_hi),
            "bar_lo_outer": float(bar_lo_outer),
            "bar_hi_outer": float(bar_hi_outer),
            "px": [float(xm), float(bar_top)],
            "px_lo": float(uc),
            "px_hi": float(lc),
            "note": sem_note,
        })
    return out


# ---------------------------------------------------------------------------
# Figure 3B (page 4, 900 dpi): Pyrene fluorescence timecourse
# ---------------------------------------------------------------------------

def calibrate_fig3b(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    # X ticks: 0 (cL) to 30 (cR), major every 5
    t_ticks = find_inward_ticks(img, box, "top", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    x_marks = [(cL, 0.0), (cR, 30.0)]
    for val in (5.0, 10.0, 15.0, 20.0, 25.0):
        ep = cL + w * val / 30.0
        near = [p for p, _ in t_ticks if abs(p - ep) <= 15]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y ticks: 0 (rB) up to 45 (rT); major ticks at 10, 20, 30, 40
    l_ticks = find_inward_ticks(img, box, "left", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    y_marks = [(rB, 0.0)]
    for val in (10.0, 20.0, 30.0, 40.0):
        ep = rB - h * val / 45.0
        near = [p for p, _ in l_ticks if abs(p - ep) <= 20]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_fig3b(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_fig3b(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "actin_alone":
        mask = (gray < 100) & (sat < 35)
    elif series == "plus_plk":
        mask = (b - r > 80) & (g - r > 50)
    elif series == "plus_coacervate":
        mask = (r - g > 80) & (r - b > 80)
    else:
        raise ValueError(series)

    # Mask legend at top (y > 35.5)
    mask[: int(round(ya.to_pixel(35.5))), :] = False
    # Exclude bottom ticks (major ticks extend ~36 px inward); keep side inset narrow so t=29.5 min is read
    if series == "actin_alone":
        mask[int(rB) - 46 :, :] = False
        mask[:, : int(cL) + 12] = False
        mask[:, int(cR) - 12 :] = False

    return _trace_mask_at(mask, xa, ya, box, at, half_w=3, fallback_w=10, inset=3, extra_px=1.2)


# ---------------------------------------------------------------------------
# Figure 3C (page 4, 900 dpi): Log-log assembly rate vs [pLK_126]
# ---------------------------------------------------------------------------

def calibrate_fig3c(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40)
    cL, rT, cR, rB = box
    # X log ticks: 0.1 (cL), 1, 10, 100, 1000
    t_ticks = find_inward_ticks(img, box, "top", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    maj_x = sorted([p for p, d in t_ticks if p > cL + 100 and p < cR - 50])
    # Expect 4 interior major decade ticks: 1, 10, 100, 1000
    if len(maj_x) < 4:
        raise RuntimeError(f"Fig 3C x major ticks: {maj_x}")
    # Keep the 4 evenly spaced decade ticks
    x_marks = [(cL, 0.1)] + [(maj_x[i], 10.0**i) for i in range(4)]
    xa = D.Axis("log10", x_marks)

    # Y log ticks: 0.1 and 1.0 major ticks, plus minor ticks (0.02..0.09, 0.2..0.9, 2..4)
    l_ticks = find_inward_ticks(img, box, "left", depth_start=6, depth_end=30, min_depth=10, max_width=14)
    maj_y = sorted([p for p, d in l_ticks if d >= 18 and rT + 50 < p < rB - 50])
    if len(maj_y) < 2:
        raise RuntimeError(f"Fig 3C y major ticks: {maj_y}")
    ya0 = D.Axis("log10", [(maj_y[0], 1.0), (maj_y[-1], 0.1)])
    y_marks = [(maj_y[0], 1.0), (maj_y[-1], 0.1)]
    for val in (0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2.0, 3.0, 4.0):
        ep = ya0.to_pixel(val)
        if rT + 15 < ep < rB - 15:
            near = [p for p, _ in l_ticks if abs(p - ep) <= 10]
            if near:
                y_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    ya = D.Axis("log10", sorted(y_marks))
    return xa, ya, box


def read_fig3c(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_fig3c(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    red = (r - g > 80) & (r - b > 80)
    cyan = (b - r > 80) & (g - r > 50)

    # Horizontal dashed control baselines (bottom region y < 0.035)
    r_base_top = int(round(ya.to_pixel(0.033)))
    if series in ("control_coacervate", "control_plk"):
        m = red if series == "control_coacervate" else cyan
        sub = m[r_base_top : int(rB) - 5, int(cL) + 15 : int(cR) - 15]
        ys, xs = np.where(sub)
        if len(ys) == 0:
            raise RuntimeError(f"Fig 3C {series} dashed line not found")
        row = float(np.median(ys)) + r_base_top
        h_px = max(2.0, 0.5 * (float(np.percentile(ys, 90)) - float(np.percentile(ys, 10))) + ya.residual_px + 1.5)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, row, h_px)
        out = []
        for xv in at:
            col = xa.to_pixel(xv)
            out.append({
                "x": float(xv),
                "y": float(y),
                "lo": float(lo),
                "hi": float(hi),
                "px": [float(col), float(row)],
                "px_lo": float(px_lo),
                "px_hi": float(px_hi),
                "note": f"horizontal dashed control, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            })
        return out

    # Mask borders, dashed baselines, and legend region (x > 40 uM, y < 0.4)
    for m in (red, cyan):
        m[: int(rT) + 6, :] = False
        m[r_base_top:, :] = False
        m[:, : int(cL) + 6] = False
        m[:, int(cR) - 6 :] = False
        m[int(round(ya.to_pixel(0.45))) :, int(round(xa.to_pixel(45.0))) :] = False

    # Detect all 13 red circle centers via binary opening + distance transform peaks
    yy, xx = np.ogrid[-10:11, -10:11]
    disk = (xx * xx + yy * yy) <= 10 * 10
    red_op = ndimage.binary_opening(red, structure=disk)
    r_dist = ndimage.distance_transform_edt(red_op)
    r_ds = ndimage.gaussian_filter(r_dist, sigma=2.0)
    r_mx = ndimage.maximum_filter(r_ds, size=29)
    ry, rx = np.where((r_ds == r_mx) & (r_dist >= 18.0))
    red_pts = sorted([(float(x), float(y), float(r_dist[y, x])) for y, x in zip(ry, rx)])

    if series == "coacervate":
        out = []
        for cx, cy, rad in red_pts:
            xv = xa.to_value(cx)
            y, lo, hi, px_lo, px_hi = _interval_from_px(ya, cy, max(2.0, ya.residual_px + 2.0))
            out.append({
                "x": float(xv),
                "y": float(y),
                "lo": float(lo),
                "hi": float(hi),
                "px": [cx, cy],
                "px_lo": float(px_lo),
                "px_hi": float(px_hi),
                "note": f"red circle r={rad:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            })
        return out

    if series == "plk_alone":
        # Cyan squares share the first 11 x positions of the red circles (0.3 .. 30 uM);
        # determine canonical square height from the unobstructed 30 uM square after morphological opening
        cyan_op = ndimage.binary_opening(cyan, structure=disk)
        sq_h = 56.0
        c30 = int(round(red_pts[10][0]))
        sub30 = cyan_op[:, max(int(cL) + 6, c30 - 45) : min(int(cR) - 6, c30 + 45)]
        ys30, xs30 = np.where(sub30)
        if len(ys30) > 1000:
            sq_h = float(ys30.max() - ys30.min() + 1)
        half_sq = 0.5 * sq_h

        out = []
        # Check first 11 concentrations (<= 40 uM)
        for rx_px, ry_px, _ in red_pts:
            if xa.to_value(rx_px) > 40.0:
                continue
            c_win0 = max(int(cL) + 6, int(round(rx_px - half_sq - 4)))
            c_win1 = min(int(cR) - 6, int(round(rx_px + half_sq + 4)))
            # Look in vertical strip around rx_px
            # Exclude thin diagonal connecting lines by looking for rows where cyan+red width >= 0.75*sq_h
            # or pure cyan width >= 0.45*sq_h (for partially occluded squares)
            strip_c = cyan[:, c_win0:c_win1]
            strip_r = red[:, c_win0:c_win1]
            # Restrict vertically to plausible region around square (within +-250 px of red point, below y=2)
            y_lo_lim = max(int(rT) + 6, int(round(ry_px)) - 120)
            y_hi_lim = min(r_base_top - 5, int(round(ry_px)) + 260)
            sc = strip_c[y_lo_lim:y_hi_lim, :]
            sr = strip_r[y_lo_lim:y_hi_lim, :]
            c_row_w = sc.sum(axis=1)
            # Find contiguous vertical band of cyan belonging to the square
            # For unobstructed squares, c_row_w has a plateau ~56 px wide over ~56 rows.
            # For squares occluded by a red circle (0.3, 0.375, 0.5, 15 uM), look at
            # columns just inside the square's left/right edges (rx_px +- 24) or top/bottom edges!
            # Check vertical extent where cyan exists across a horizontal width >= 12 px
            strong_rows = np.where(c_row_w >= 14)[0] + y_lo_lim
            if len(strong_rows) == 0:
                continue
            # Also check if connecting line extends vertically: trim thin tails where c_row_w < 24
            wide_rows = np.where(c_row_w >= 24)[0] + y_lo_lim
            # Wait: at 15 uM where red circle covers the middle of the cyan square, the cyan corners
            # at top and bottom still have width ~56 px, while the middle rows have only thin left/right slivers (~8 px)!
            # So look at BOTH top and bottom cyan horizontal edge spans!
            # Narrow window specifically to the square's expected x-span [rx_px - half_sq + 2 .. rx_px + half_sq - 2]
            x_in0 = max(int(cL) + 6, int(round(rx_px - half_sq + 2)))
            x_in1 = min(int(cR) - 6, int(round(rx_px + half_sq - 2)))
            in_c = cyan[y_lo_lim:y_hi_lim, x_in0:x_in1]
            in_r = red[y_lo_lim:y_hi_lim, x_in0:x_in1]
            # A row belongs to the cyan square if cyan pixels span wide (>= 26 px) OR cyan appears on both left & right edges
            row_w = in_c.sum(axis=1)
            both_edges = in_c[:, :10].any(axis=1) & in_c[:, -10:].any(axis=1)
            sq_rows = np.where((row_w >= 24) | both_edges)[0] + y_lo_lim
            if len(sq_rows) == 0:
                continue
            # Keep the cluster closest to the square's vertical center (if multiple due to touching neighbor)
            # Note: neighbors at 0.3 / 0.375 / 0.5 overlap slightly in x! Use narrow central columns [rx_px-12 : rx_px+12]
            cx0 = max(int(cL) + 6, int(round(rx_px - 12)))
            cx1 = min(int(cR) - 6, int(round(rx_px + 12)))
            mid_c = cyan[y_lo_lim:y_hi_lim, cx0:cx1]
            mid_r = red[y_lo_lim:y_hi_lim, cx0:cx1]
            mid_w = mid_c.sum(axis=1)
            # Clean rows in central strip where cyan covers >= 12 of 24 cols (not thin 10px connecting line)
            # Wait: thin connecting line is ~9 px thick vertically at diagonal -> covers < 15 rows; square covers ~56 rows!
            # Better: look at cyan pixels at columns [rx_px-15 .. rx_px+15] and find the top and bottom edges of the square!
            c_rows_mid = np.where(mid_w >= 12)[0] + y_lo_lim
            if len(c_rows_mid) == 0:
                # Completely covered in middle (like 15 uM): use full width edges!
                y_top = float(sq_rows.min())
                y_bot = float(sq_rows.max())
                cy = 0.5 * (y_top + y_bot)
                occ = "4-corner bbox (red overlap)"
            else:
                # Largest cluster in c_rows_mid
                cls = _cluster_1d(c_rows_mid, gap=4)
                # At 15 uM there are two thin top/bottom strips separated by red circle!
                if len(cls) >= 2 and (cls[-1][1] - cls[0][0]) <= sq_h + 6:
                    y_top, y_bot = float(cls[0][0]), float(cls[-1][1])
                    cy = 0.5 * (y_top + y_bot)
                    occ = "top+bottom edges (red center overlap)"
                else:
                    # Pick cluster with largest vertical height (not a diagonal line)
                    # Wait: check if red sits immediately below y_bot or above y_top
                    best_c = max(cls, key=lambda g: g[1] - g[0])
                    y_top, y_bot = float(best_c[0]), float(best_c[1])
                    span = y_bot - y_top + 1
                    if span >= 0.85 * sq_h:
                        cy = 0.5 * (y_top + y_bot)
                        occ = f"unobstructed (h={span:.0f}px)"
                    else:
                        # Check which side is occluded by red!
                        red_above = red[max(0, int(y_top) - 10) : int(y_top) + 2, cx0:cx1].sum()
                        red_below = red[int(y_bot) - 1 : min(img.shape[0], int(y_bot) + 11), cx0:cx1].sum()
                        if red_below > red_above:
                            cy = y_top + half_sq - 0.5
                            occ = f"top-edge anchored (bottom red-occluded, vis={span:.0f}px)"
                        else:
                            cy = y_bot - half_sq + 0.5
                            occ = f"bottom-edge anchored (top red-occluded, vis={span:.0f}px)"

            xv = xa.to_value(rx_px)
            err_px = 3.0 if "occluded" in occ or "overlap" in occ else 2.0
            y, lo, hi, px_lo, px_hi = _interval_from_px(ya, cy, max(err_px, ya.residual_px + err_px))
            out.append({
                "x": float(xv),
                "y": float(y),
                "lo": float(lo),
                "hi": float(hi),
                "px": [float(rx_px), float(cy)],
                "px_lo": float(px_lo),
                "px_hi": float(px_hi),
                "note": f"cyan square ({occ}), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            })
        return out

    raise ValueError(series)


# ---------------------------------------------------------------------------
# Figure 3D (page 4, 900 dpi): Pyrene excimer fluorescence timecourse
# ---------------------------------------------------------------------------

def calibrate_fig3d(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40, min_h_frac=0.55, min_v_frac=0.55)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    # X ticks: 0 (cL) to 30 (cR), major every 5
    t_ticks = find_inward_ticks(img, box, "top", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    x_marks = [(cL, 0.0), (cR, 30.0)]
    for val in (5.0, 10.0, 15.0, 20.0, 25.0):
        ep = cL + w * val / 30.0
        near = [p for p, _ in t_ticks if abs(p - ep) <= 15]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y ticks: 10 (rB) to 37.5 (rT); major right ticks at 15, 20, 25, 30, 35
    r_ticks = find_inward_ticks(img, box, "right", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    y_marks = [(rB, 10.0)]
    for val in (15.0, 20.0, 25.0, 30.0, 35.0):
        ep = rB - h * (val - 10.0) / 27.5
        near = [p for p, _ in r_ticks if abs(p - ep) <= 20]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_fig3d(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_fig3d(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "actin_alone":
        mask = (gray < 100) & (sat < 35)
    elif series == "plus_plk":
        mask = (b - r > 80) & (g - r > 50)
    elif series == "plus_coacervate":
        mask = (r - g > 80) & (r - b > 80)
    else:
        raise ValueError(series)

    # Mask legend at top-center (x > 2.5 min, y > 28)
    mask[: int(round(ya.to_pixel(28.0))), int(round(xa.to_pixel(2.5))) :] = False
    if series == "actin_alone":
        mask[int(rB) - 46 :, :] = False
        mask[: int(rT) + 46, :] = False
        mask[:, : int(cL) + 12] = False
        mask[:, int(cR) - 12 :] = False

    return _trace_mask_at(mask, xa, ya, box, at, half_w=3, fallback_w=10, inset=3, extra_px=1.2)


# ---------------------------------------------------------------------------
# Figure 4D (page 5, 1200 dpi): Semi-log TMR-actin intensity line scans
# ---------------------------------------------------------------------------

def calibrate_fig4d(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40, min_h_frac=0.60, min_v_frac=0.55)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    # X ticks: -4 (cL) to +4 (cR), major ticks every 1 um
    t_ticks = find_inward_ticks(img, box, "top", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    x_marks = [(cL, -4.0), (cR, 4.0)]
    for val in (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0):
        ep = cL + w * (val + 4.0) / 8.0
        near = [p for p, _ in t_ticks if abs(p - ep) <= 18]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y log ticks: 10^0 (rB), 10^1, 10^2, 10^3, 10^4 (major right ticks)
    r_ticks = find_inward_ticks(img, box, "right", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    maj_y = sorted([p for p, d in r_ticks if rT + 30 < p < rB - 30])
    # Expect 4 interior major ticks: 10^4, 10^3, 10^2, 10^1
    if len(maj_y) < 4:
        raise RuntimeError(f"Fig 4D y major ticks: {maj_y}")
    y_marks = [(maj_y[0], 1e4), (maj_y[1], 1e3), (maj_y[2], 1e2), (maj_y[3], 1e1), (rB, 1.0)]
    ya = D.Axis("log10", sorted(y_marks))
    return xa, ya, box


_FIG4D_COLORS = {
    "1p5uM": ((144, 39, 142), 62.0),   # purple
    "0p5uM": ((15, 172, 237), 65.0),   # blue/cyan
    "0p25uM": ((15, 164, 82), 62.0),   # green
    "0p10uM": ((254, 240, 2), 68.0),   # yellow
    "0p05uM": ((245, 145, 31), 58.0),  # orange
    "0p01uM": ((235, 29, 35), 62.0),   # red
}


def read_fig4d(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_fig4d(img)
    cL, rT, cR, rB = box
    rgb, tol = _FIG4D_COLORS[series]
    mask = _rgb_dist_mask(img, rgb, tol)

    # Mask out legend box at bottom-center (x in [-3.75, 1.9], y < 6.5)
    lc0 = int(round(xa.to_pixel(-3.78)))
    lc1 = int(round(xa.to_pixel(1.92)))
    lr0 = int(round(ya.to_pixel(6.5)))
    mask[lr0:, lc0:lc1] = False

    # Use slightly wider window (open circles are ~20 px wide) so we find center of circle/line cleanly
    return _trace_mask_at(mask, xa, ya, box, at, half_w=5, fallback_w=16, inset=3, extra_px=2.0, mode="mean")


# ---------------------------------------------------------------------------
# Figure 4E (page 5, 1200 dpi): Log-log <Intensity>_droplet vs [TMR-actin]
# ---------------------------------------------------------------------------

def calibrate_fig4e(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40, min_h_frac=0.60, min_v_frac=0.55)
    cL, rT, cR, rB = box
    # X log major ticks on top spine: 0.01, 0.1, 1.0
    t_ticks = find_inward_ticks(img, box, "top", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    maj_x = sorted([p for p, d in t_ticks if cL + 30 < p < cR - 30])
    if len(maj_x) < 3:
        raise RuntimeError(f"Fig 4E x major ticks: {maj_x}")
    xa = D.Axis("log10", [(maj_x[0], 0.01), (maj_x[1], 0.1), (maj_x[2], 1.0)])

    # Y log major ticks on right spine: 10^4, 10^3, 10^2, 10^1 (rB)
    r_ticks = find_inward_ticks(img, box, "right", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    maj_y = sorted([p for p, d in r_ticks if rT + 30 < p < rB - 30])
    if len(maj_y) < 3:
        raise RuntimeError(f"Fig 4E y major ticks: {maj_y}")
    ya = D.Axis("log10", [(maj_y[0], 1e4), (maj_y[1], 1e3), (maj_y[2], 1e2), (rB, 1e1)])
    return xa, ya, box


def read_fig4e(img: np.ndarray, series: str) -> list[dict]:
    xa, ya, box = calibrate_fig4e(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "periphery":
        mask = (b - r > 80) & (g - r > 50)
    elif series == "interior":
        mask = (gray < 85) & (sat < 30)
    elif series == "exterior":
        mask = (r - g > 80) & (r - b > 80)
    else:
        raise ValueError(series)

    # Mask frame borders (tight margin for colored markers so 0.01 uM red square near y=10 isn't clipped)
    m_in = 14 if series == "interior" else 4
    mask[: int(rT) + m_in, :] = False
    mask[int(rB) - m_in :, :] = False
    mask[:, : int(cL) + m_in] = False
    mask[:, int(cR) - m_in :] = False

    # Mask legend box at top-left (x in [0.011, 0.075], y in [1400, 14000])
    lc0 = int(round(xa.to_pixel(0.011)))
    lc1 = int(round(xa.to_pixel(0.078)))
    lr0 = int(round(ya.to_pixel(15000.0)))
    lr1 = int(round(ya.to_pixel(1300.0)))
    mask[lr0:lr1, lc0:lc1] = False

    # Detect marker at each of the 6 experimental concentrations
    xs_exp = [0.01, 0.05, 0.1, 0.25, 0.5, 1.5]
    out = []
    for xv in xs_exp:
        cx0 = int(round(xa.to_pixel(xv)))
        a = max(int(cL) + m_in, cx0 - 45)
        b_col = min(int(cR) - m_in, cx0 + 45)
        sub = mask[:, a:b_col]
        lab, n = ndimage.label(sub)
        best = None
        for k in range(1, n + 1):
            ys, xs = np.where(lab == k)
            if len(xs) >= 120:
                if best is None or len(xs) > best[0]:
                    best = (len(xs), ys, xs + a)
        if best is None:
            out.append({"x": float(xv), "y": None, "found": False, "note": "marker not found"})
            continue
        _, ys, xs = best
        cx = float(xs.mean())
        # For squares/circles mean == bbox center; for triangles, mean is centroid (matches 4F triangles too)
        cy = float(ys.mean())
        x_read = xa.to_value(cx)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, cy, max(2.5, ya.residual_px + 2.5))
        out.append({
            "x": float(xv),  # nominal experimental [TMR-actin] (verified x_read within 1%)
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "px": [cx, cy],
            "px_lo": float(px_lo),
            "px_hi": float(px_hi),
            "note": f"x_img={x_read:.4g} uM, area={len(xs)}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


# ---------------------------------------------------------------------------
# Figure 4F (page 5, 1200 dpi): Log-log Auxiliary RAM & Filament Boundary
# ---------------------------------------------------------------------------

def calibrate_fig4f(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40, min_h_frac=0.60, min_v_frac=0.55)
    cL, rT, cR, rB = box
    t_ticks = find_inward_ticks(img, box, "top", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    maj_x = sorted([p for p, d in t_ticks if cL + 30 < p < cR - 30])
    if len(maj_x) < 3:
        raise RuntimeError(f"Fig 4F x major ticks: {maj_x}")
    xa = D.Axis("log10", [(maj_x[0], 0.01), (maj_x[1], 0.1), (maj_x[2], 1.0)])

    r_ticks = find_inward_ticks(img, box, "right", depth_start=6, depth_end=30, min_depth=10, max_width=14)
    maj_y = sorted([p for p, d in r_ticks if d >= 18 and rT + 50 < p < rB - 50])
    if len(maj_y) < 1:
        raise RuntimeError(f"Fig 4F y major ticks: {r_ticks}")
    ya0 = D.Axis("log10", [(rT, 100.0), (maj_y[0], 10.0), (rB, 1.0)])
    y_marks = [(rT, 100.0), (maj_y[0], 10.0), (rB, 1.0)]
    for val in (2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0):
        ep = ya0.to_pixel(val)
        near = [p for p, _ in r_ticks if abs(p - ep) <= 10]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    ya = D.Axis("log10", sorted(y_marks))
    return xa, ya, box


def read_fig4f(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_fig4f(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    # Legend box at bottom-right: x > 0.11 uM, y < 3.8
    lc0 = int(round(xa.to_pixel(0.115)))
    lr0 = int(round(ya.to_pixel(3.8)))

    if series == "filament_boundary":
        red = (r - g > 80) & (r - b > 80)
        red[: int(rT) + 6, :] = False
        red[int(rB) - 6 :, :] = False
        red[:, : int(cL) + 6] = False
        red[:, int(cR) - 6 :] = False
        red[lr0:, lc0:] = False
        # Fit straight line in pixel space to all detected red dashes of the Filament Boundary
        ys, xs = np.where(red)
        if len(xs) < 30:
            raise RuntimeError("Fig 4F red filament boundary line not found")
        slope, intercept = np.polyfit(xs.astype(float), ys.astype(float), 1)
        resid_line = float(np.percentile(np.abs(ys - (slope * xs + intercept)), 90))
        out = []
        for xv in at:
            col = xa.to_pixel(xv)
            row = slope * col + intercept
            # Also verify local red dash density near this (col, row)
            h_px = max(2.5, resid_line + ya.residual_px + 1.5)
            y, lo, hi, px_lo, px_hi = _interval_from_px(ya, row, h_px)
            out.append({
                "x": float(xv),
                "y": float(y),
                "lo": float(lo),
                "hi": float(hi),
                "px": [float(col), float(row)],
                "px_lo": float(px_lo),
                "px_hi": float(px_hi),
                "note": f"red dashed line fit (p90_res={resid_line:.2f}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            })
        return out

    if series == "ram_periphery":
        mask = (b - r > 80) & (g - r > 50)
    elif series == "ram_interior":
        mask = (gray < 85) & (sat < 30)
    else:
        raise ValueError(series)

    m_in = 14 if series == "ram_interior" else 6
    mask[: int(rT) + m_in, :] = False
    mask[int(rB) - m_in :, :] = False
    mask[:, : int(cL) + m_in] = False
    mask[:, int(cR) - m_in :] = False
    mask[lr0:, lc0:] = False

    xs_exp = [0.01, 0.05, 0.1, 0.25, 0.5, 1.5]
    out = []
    for xv in xs_exp:
        cx0 = int(round(xa.to_pixel(xv)))
        a = max(int(cL) + m_in, cx0 - 45)
        b_col = min(int(cR) - m_in, cx0 + 45)
        sub = mask[:, a:b_col]
        lab, n = ndimage.label(sub)
        best = None
        for k in range(1, n + 1):
            ys, xs = np.where(lab == k)
            if len(xs) >= 120:
                if best is None or len(xs) > best[0]:
                    best = (len(xs), ys, xs + a)
        if best is None:
            out.append({"x": float(xv), "y": None, "found": False, "note": "marker not found"})
            continue
        _, ys, xs = best
        cx = float(xs.mean())
        cy = float(ys.mean())
        x_read = xa.to_value(cx)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, cy, max(2.5, ya.residual_px + 2.5))
        out.append({
            "x": float(xv),
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "px": [cx, cy],
            "px_lo": float(px_lo),
            "px_hi": float(px_hi),
            "note": f"x_img={x_read:.4g} uM, area={len(xs)}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


# ---------------------------------------------------------------------------
# Figure S1B (page 25, 600 dpi): Aspect ratio during droplet coalescence
# ---------------------------------------------------------------------------

def calibrate_figS1b(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40)
    cL, rT, cR, rB = box
    h = rB - rT
    # Top ticks: 7 evenly spaced major ticks at 0.00 .. 0.06 s
    t_ticks = find_inward_ticks(img, box, "top", depth_start=6, depth_end=32, min_depth=16, max_width=14)
    xs_t = sorted([p for p, _ in t_ticks])
    if len(xs_t) != 7:
        raise RuntimeError(f"Fig S1B expected 7 x-ticks, found {xs_t}")
    xa = D.Axis("linear", [(xs_t[i], i * 0.01) for i in range(7)])

    # Left major ticks at 1.8 (rT), 1.6, 1.4, 1.2, 1.0, 0.8 (rB)
    l_ticks = find_inward_ticks(img, box, "left", depth_start=6, depth_end=32, min_depth=16, max_width=14)
    y_marks = [(rT, 1.8), (rB, 0.8)]
    for val in (1.6, 1.4, 1.2, 1.0):
        ep = rT + h * (1.8 - val) / 1.0
        near = [p for p, _ in l_ticks if abs(p - ep) <= 18]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_figS1b(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_figS1b(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "fit":
        red = (r - g > 60) & (r - b > 60)
        # Mask legend at top-right (t > 0.032 s, y > 1.45)
        red[: int(round(ya.to_pixel(1.45))), int(round(xa.to_pixel(0.032))) :] = False
        return _trace_mask_at(red, xa, ya, box, at, half_w=3, fallback_w=10, inset=4, extra_px=1.8)

    if series == "data":
        dark = (gray < 110) & (sat < 45)
        red_curve = (r - g > 60) & (r - b > 60)
        for m in (dark, red_curve):
            m[: int(rT) + 30, :] = False
            m[int(rB) - 30 :, :] = False
            m[:, : int(cL) + 30] = False
            m[:, int(cR) - 30 :] = False
            # Mask legend
            m[: int(round(ya.to_pixel(1.45))), int(round(xa.to_pixel(0.032))) :] = False

        out = []
        for tv in [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]:
            cx0 = int(round(xa.to_pixel(tv)))
            # Search vertical strip +-35 px around tick column
            w0, w1 = cx0 - 35, cx0 + 36
            sub = dark[:, w0:w1]
            row_w = sub.sum(axis=1)
            # Find vertical extent of dark error bar stem (allow red fit line crossing to bridge stem)
            stem = (dark | red_curve)[:, cx0 - 4 : cx0 + 5].any(axis=1)
            stem_rows = np.where(stem)[0]
            if len(stem_rows) == 0:
                out.append({"x": tv, "y": None, "found": False, "note": "not found"})
                continue
            cls = _cluster_1d(stem_rows, gap=28)
            s_top, s_bot, _ = max(cls, key=lambda g: g[1] - g[0])
            # Trim any non-dark red-only rows at extreme ends of [s_top, s_bot]
            dark_stem = np.where(dark[s_top : s_bot + 1, cx0 - 12 : cx0 + 13].any(axis=1))[0] + s_top
            s_top, s_bot = int(dark_stem[0]), int(dark_stem[-1])

            # Marker center: thickest solid square block (rows where narrow central strip cx0-9..cx0+9 has >=14 dark px)
            core = dark[s_top : s_bot + 1, cx0 - 9 : cx0 + 10].sum(axis=1)
            # Also find the dx bar (widest rows in [s_top..s_bot])
            mid_w = row_w[s_top : s_bot + 1]
            marker_rows = np.where(core >= 14)[0] + s_top
            if len(marker_rows) == 0:
                marker_rows = np.where(mid_w >= 0.7 * mid_w.max())[0] + s_top
            # Exclude top/bottom caps if they're >= 14 px wide by restricting to inner 60% of [s_top, s_bot]
            inner = marker_rows[(marker_rows > s_top + 25) & (marker_rows < s_bot - 25)]
            cy = float(np.median(inner)) if len(inner) else float(np.median(marker_rows))

            # Top cap center: dark cap rows within [s_top .. s_top + 12]
            top_cap_rows = np.where(row_w[s_top : s_top + 12] >= 12)[0] + s_top
            bot_cap_rows = np.where(row_w[max(0, s_bot - 11) : s_bot + 1] >= 12)[0] + max(0, s_bot - 11)
            uc = float(np.mean(top_cap_rows)) if len(top_cap_rows) else float(s_top)
            lc = float(np.mean(bot_cap_rows)) if len(bot_cap_rows) else float(s_bot)

            y, lo, hi, px_lo, px_hi = _interval_from_px(ya, cy, max(3.0, ya.residual_px + 2.5))
            bar_hi = ya.to_value(uc)
            bar_lo = ya.to_value(lc)
            end_err = max(2.5, ya.residual_px + 2.0)
            bar_hi_outer = ya.to_value(uc - end_err)
            bar_lo_outer = ya.to_value(lc + end_err)
            out.append({
                "x": float(tv),
                "y": float(y),
                "lo": float(lo),
                "hi": float(hi),
                "bar_lo": float(min(bar_lo, y)),
                "bar_hi": float(max(bar_hi, y)),
                "bar_lo_outer": float(min(bar_lo_outer, bar_lo, y)),
                "bar_hi_outer": float(max(bar_hi_outer, bar_hi, y)),
                "px": [float(cx0), float(cy)],
                "px_lo": float(uc),
                "px_hi": float(lc),
                "note": f"dy caps=({uc:.1f},{lc:.1f}), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            })
        return out

    raise ValueError(series)


# ---------------------------------------------------------------------------
# Figure S2C (page 26, 600 dpi): Stacked Actin+DIC (top) and BSA+DIC (bottom)
# ---------------------------------------------------------------------------

def calibrate_figS2c(img: np.ndarray, panel: str):
    gray, sat = _gray_sat(img)
    dark = (gray < 175) & (sat < 40)
    H, W = dark.shape
    # Find 3 horizontal spines (top, middle divider, bottom) via vectorized mean
    row_frac = dark.mean(axis=1)
    hc = _cluster_1d(np.where(row_frac >= 0.45)[0], gap=8)
    h_scores = [float(row_frac[a : b + 1].max()) for a, b, _ in hc]
    max_h = max(h_scores)
    h_sp = [float(c[2]) for c, sc in zip(hc, h_scores) if sc >= 0.78 * max_h]
    if len(h_sp) < 3:
        raise RuntimeError(f"Fig S2C expected 3 H-spines, found {h_sp}")
    rTop, rMid, rBot = h_sp[0], h_sp[len(h_sp) // 2], h_sp[-1]

    # Left & right vertical spines spanning rTop..rBot
    ra, rb = int(round(rTop)), int(round(rBot))
    col_frac = dark[ra:rb].mean(axis=0)
    vc = _cluster_1d(np.where(col_frac >= 0.55)[0], gap=8)
    cL, cR = float(vc[0][2]), float(vc[-1][2])

    rT, rB = (rTop, rMid) if panel == "top" else (rMid, rBot)
    box = (cL, rT, cR, rB)

    # X axis: -4 (cL) to +4 (cR), check top ticks at -2, 0, +2
    w = cR - cL
    ticks_x = find_inward_ticks(img, (cL, rTop, cR, rBot), "top", dark_thr=175, depth_start=3, depth_end=20, min_depth=5, max_width=10)
    x_marks = [(cL, -4.0), (cR, 4.0)]
    for val in (-2.0, 0.0, 2.0):
        ep = cL + w * (val + 4.0) / 8.0
        near = [p for p, _ in ticks_x if abs(p - ep) <= 12]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y axis: inward left/right ticks for 1.0, 0.5, 0.0 (slightly inset from rT and rB)
    l_ticks = find_inward_ticks(img, box, "left", dark_thr=175, depth_start=3, depth_end=20, min_depth=6, max_width=10, margin=2)
    r_ticks = find_inward_ticks(img, box, "right", dark_thr=175, depth_start=3, depth_end=20, min_depth=6, max_width=10, margin=2)
    h = rB - rT
    y_marks = []
    for val, frac in ((1.0, 0.016), (0.5, 0.50), (0.0, 0.98)):
        ep = rT + frac * h
        cands = [p for p, _ in l_ticks + r_ticks if abs(p - ep) <= 14]
        if cands:
            y_marks.append((float(np.median(cands)), val))
    if len(y_marks) < 2:
        y_marks = [(rT + 0.016 * h, 1.0), (rB - 0.02 * h, 0.0)]
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_figS2c(img: np.ndarray, at: list[float], panel: str, series: str) -> list[dict]:
    xa, ya, box = calibrate_figS2c(img, panel)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "actin":
        mask = (g - r > 35) & (g - b > 15) & (g > 80)
    elif series == "bsa":
        mask = (r - b > 70) & (r > 160) & (g > 90)
    elif series == "dic":
        mask = (gray < 110) & (sat < 35)
    else:
        raise ValueError(series)

    # Mask legend box inside each subpanel:
    # Top legend: x in [-1.2, 1.3], y in [0.64, 0.96]; Bottom legend: x in [-1.2, 1.3], y in [0.52, 0.80]
    lc0 = int(round(xa.to_pixel(-1.25)))
    lc1 = int(round(xa.to_pixel(1.35)))
    if panel == "top":
        lr0 = int(round(ya.to_pixel(0.96)))
        lr1 = int(round(ya.to_pixel(0.62)))
    else:
        lr0 = int(round(ya.to_pixel(0.82)))
        lr1 = int(round(ya.to_pixel(0.50)))
    mask[lr0:lr1, lc0:lc1] = False

    if series == "dic":
        mask[int(rB) - 18 :, :] = False
        mask[: int(rT) + 18, :] = False
        mask[:, : int(cL) + 18] = False
        mask[:, int(cR) - 18 :] = False

    return _trace_mask_at(mask, xa, ya, box, at, half_w=3, fallback_w=10, inset=2, extra_px=2.0)


# ---------------------------------------------------------------------------
# Figure S3B (page 27, 600 dpi): 647-BSA flow-in RAM vs Time
# ---------------------------------------------------------------------------

def calibrate_figS3b(img: np.ndarray):
    box = find_frame(img, dark_thr=205, sat_thr=35)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    # X ticks: -10 (cL) to 70 (cR) every 10 s
    b_ticks = find_inward_ticks(img, box, "bottom", dark_thr=180, depth_start=4, depth_end=24, min_depth=6, max_width=10)
    x_marks = [(cL, -10.0), (cR, 70.0)]
    for val in (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0):
        ep = cL + w * (val + 10.0) / 80.0
        near = [p for p, _ in b_ticks if abs(p - ep) <= 14]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y ticks: 0 (rB) to 10 (rT) every 1 unit
    r_ticks = find_inward_ticks(img, box, "right", dark_thr=180, depth_start=4, depth_end=24, min_depth=6, max_width=10)
    y_marks = [(rT, 10.0), (rB, 0.0)]
    for val in range(1, 10):
        ep = rB - h * val / 10.0
        near = [p for p, _ in r_ticks if abs(p - ep) <= 14]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), float(val)))
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_figS3b(img: np.ndarray) -> list[dict]:
    xa, ya, box = calibrate_figS3b(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray = (r + g + b) / 3.0

    # Dots are black, green (t=0), or red (t=35) filled circles with dark outline
    dot_mask = (gray < 105) | ((g - r > 40) & (g - b > 20)) | ((r - g > 50) & (r - b > 50))
    dot_mask[: int(rT) + 20, :] = False
    dot_mask[int(rB) - 20 :, :] = False
    dot_mask[:, : int(cL) + 15] = False
    dot_mask[:, int(cR) - 15 :] = False
    # Exclude "Flow-in BSA" label in blue shaded band (t in [2, 33], y > 7.5)
    dot_mask[: int(round(ya.to_pixel(7.2))), int(round(xa.to_pixel(2.0))) : int(round(xa.to_pixel(33.5)))] = False

    # Use distance transform + local maxima to cleanly split touching dots (e.g. t=-8 and -6 s)
    dist = ndimage.distance_transform_edt(dot_mask)
    dist_s = ndimage.gaussian_filter(dist, sigma=2.0)
    mx = ndimage.maximum_filter(dist_s, size=25)
    py, px = np.where((dist_s == mx) & (dist >= 10.0))

    pts = sorted([(float(x), float(y), float(dist[y, x])) for y, x in zip(py, px)])
    out = []
    for cx, cy, rad in pts:
        xv = xa.to_value(cx)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, cy, max(2.5, ya.residual_px + 2.0))
        out.append({
            "x": float(xv),
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "px": [cx, cy],
            "px_lo": float(px_lo),
            "px_hi": float(px_hi),
            "note": f"marker r={rad:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


# ---------------------------------------------------------------------------
# Figure S3D (page 27, 600 dpi): TMR-actin flow-in intensity profile & baseline
# ---------------------------------------------------------------------------

def calibrate_figS3d(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    b_ticks = find_inward_ticks(img, box, "bottom", depth_start=6, depth_end=28, min_depth=10, max_width=12)
    x_marks = [(cL, -15.0), (cR, 15.0)]
    for val in (-10.0, -5.0, 0.0, 5.0, 10.0):
        ep = cL + w * (val + 15.0) / 30.0
        near = [p for p, _ in b_ticks if abs(p - ep) <= 14]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    r_ticks = find_inward_ticks(img, box, "right", depth_start=6, depth_end=28, min_depth=12, max_width=12)
    y_marks = [(rT, 700.0), (rB, 0.0)]
    for val in (600.0, 500.0, 400.0, 300.0, 200.0, 100.0):
        ep = rB - h * val / 700.0
        near = [p for p, _ in r_ticks if abs(p - ep) <= 16]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_figS3d(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_figS3d(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "baseline":
        red = (r - g > 45) & (r - b > 45)
        sub = red[int(rT) + 20 : int(rB) - 20, int(cL) + 20 : int(cR) - 20]
        ys, xs = np.where(sub)
        if len(ys) == 0:
            raise RuntimeError("Fig S3D red dashed baseline not found")
        row = float(np.median(ys)) + int(rT) + 20
        h_px = max(2.5, 0.5 * (float(np.percentile(ys, 90)) - float(np.percentile(ys, 10))) + ya.residual_px + 1.5)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, row, h_px)
        return [{
            "x": float(xv),
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(xa.to_pixel(xv)), float(row)],
            "px_lo": float(px_lo),
            "px_hi": float(px_hi),
            "note": f"red dashed exterior baseline, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        } for xv in at]

    if series == "profile":
        mask = (gray < 105) & (sat < 40)
        mask[: int(rT) + 22, :] = False
        mask[int(rB) - 22 :, :] = False
        mask[:, : int(cL) + 22] = False
        mask[:, int(cR) - 22 :] = False
        return _trace_mask_at(mask, xa, ya, box, at, half_w=3, fallback_w=10, inset=4, extra_px=2.0)

    raise ValueError(series)


# ---------------------------------------------------------------------------
# Figure S4C (page 29, 600 dpi): Post-phase-separation Actin & DIC line scans
# ---------------------------------------------------------------------------

def calibrate_figS4c(img: np.ndarray):
    box = find_frame(img, dark_thr=205, sat_thr=35)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    b_ticks = find_inward_ticks(img, box, "bottom", dark_thr=180, depth_start=4, depth_end=22, min_depth=6, max_width=10)
    x_marks = [(cL, -10.0), (cR, 10.0)]
    for val in (-5.0, 0.0, 5.0):
        ep = cL + w * (val + 10.0) / 20.0
        near = [p for p, _ in b_ticks if abs(p - ep) <= 14]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y ticks: 1.0 (~rT+14), 0.5 (~mid), 0.0 (~rB-12)
    l_ticks = find_inward_ticks(img, box, "left", dark_thr=180, depth_start=4, depth_end=24, min_depth=10, max_width=10, margin=3)
    r_ticks = find_inward_ticks(img, box, "right", dark_thr=180, depth_start=4, depth_end=24, min_depth=10, max_width=10, margin=3)
    y_marks = []
    for val, frac in ((1.0, 0.02), (0.5, 0.50), (0.0, 0.98)):
        ep = rT + frac * h
        cands = [p for p, _ in l_ticks + r_ticks if abs(p - ep) <= 18]
        if cands:
            y_marks.append((float(np.median(cands)), val))
    if len(y_marks) < 2:
        y_marks = [(rT + 0.02 * h, 1.0), (rB - 0.02 * h, 0.0)]
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_figS4c(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_figS4c(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)
    gray, sat = _gray_sat(img)

    if series == "actin":
        mask = (g - r > 30) & (g - b > 15) & (g > 75)
    elif series == "dic":
        mask = (gray < 105) & (sat < 35)
    else:
        raise ValueError(series)

    # Mask legend box at top-center (x in [-2.8, 2.8], y in [0.74, 0.97])
    lc0 = int(round(xa.to_pixel(-2.9)))
    lc1 = int(round(xa.to_pixel(2.9)))
    lr0 = int(round(ya.to_pixel(0.97)))
    lr1 = int(round(ya.to_pixel(0.72)))
    mask[lr0:lr1, lc0:lc1] = False

    if series == "dic":
        mask[int(rB) - 18 :, :] = False
        mask[: int(rT) + 18, :] = False
        mask[:, : int(cL) + 18] = False
        mask[:, int(cR) - 18 :] = False

    return _trace_mask_at(mask, xa, ya, box, at, half_w=3, fallback_w=10, inset=2, extra_px=2.0)


# ---------------------------------------------------------------------------
# Figure S4D (page 29, 600 dpi): Scatter plot RAM vs Droplet Radius (N=117)
# ---------------------------------------------------------------------------

def calibrate_figS4d(img: np.ndarray):
    box = find_frame(img, dark_thr=205, sat_thr=35, min_h_frac=0.50, min_v_frac=0.50)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    # X ticks: 1 (cL) to 5 (cR), check top ticks at 2, 3, 4
    t_ticks = find_inward_ticks(img, box, "top", dark_thr=185, depth_start=3, depth_end=16, min_depth=5, max_width=8)
    x_marks = [(cL, 1.0), (cR, 5.0)]
    for val in (2.0, 3.0, 4.0):
        ep = cL + w * (val - 1.0) / 4.0
        near = [p for p, _ in t_ticks if abs(p - ep) <= 12]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y ticks: 40 (rB) to 130 (rT)
    r_ticks = find_inward_ticks(img, box, "right", dark_thr=185, depth_start=3, depth_end=16, min_depth=4, max_width=8)
    y_marks = [(rT, 130.0), (rB, 40.0)]
    for val in range(50, 130, 10):
        ep = rB - h * (val - 40.0) / 90.0
        near = [p for p, _ in r_ticks if abs(p - ep) <= 10]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), float(val)))
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_figS4d(img: np.ndarray) -> list[dict]:
    xa, ya, box = calibrate_figS4d(img)
    cL, rT, cR, rB = box
    gray, sat = _gray_sat(img)

    # Dark circular scatter dots inside plot box
    r0, r1 = int(math.ceil(rT + 5)), int(math.floor(rB - 4))
    c0, c1 = int(math.ceil(cL + 5)), int(math.floor(cR - 5))
    sub = gray[r0:r1, c0:c1]
    mask = sub < 135
    # Exclude "r = -0.427" annotation in top-right (R > 3.55 um, RAM > 113)
    tr_c = int(round(xa.to_pixel(3.55))) - c0
    tr_r = int(round(ya.to_pixel(113.0))) - r0
    mask[: max(0, tr_r), max(0, tr_c) :] = False

    # Distance transform + local maxima separates touching dots without hallucinating false positives
    dist = ndimage.distance_transform_edt(mask)
    dist_s = ndimage.gaussian_filter(dist, sigma=1.1)
    mx = ndimage.maximum_filter(dist_s, size=9)
    py, px = np.where((dist_s == mx) & (dist >= 2.1))

    pts = sorted([(float(x + c0), float(y + r0), float(dist[y, x])) for y, x in zip(py, px)])
    out = []
    for cx, cy, rad in pts:
        xv = xa.to_value(cx)
        y, lo, hi, px_lo, px_hi = _interval_from_px(ya, cy, max(2.5, ya.residual_px + 2.0))
        out.append({
            "x": float(xv),
            "y": float(y),
            "lo": float(lo),
            "hi": float(hi),
            "px": [cx, cy],
            "px_lo": float(px_lo),
            "px_hi": float(px_hi),
            "note": f"scatter dot r={rad:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    return out


# ---------------------------------------------------------------------------
# Figure S5E (page 30, 450 dpi): Normalized FRAP of TMR-actin & 647-BSA
# ---------------------------------------------------------------------------

def calibrate_figS5e(img: np.ndarray):
    box = find_frame(img, dark_thr=150, sat_thr=40)
    cL, rT, cR, rB = box
    # Top major ticks: 0, 30, 60, 90, 120, 150, 180, 210, 240 s
    t_ticks = find_inward_ticks(img, box, "top", depth_start=5, depth_end=28, min_depth=16, max_width=10)
    maj_x = sorted([p for p, d in t_ticks if cL + 30 < p < cR - 20])
    if len(maj_x) != 9:
        raise RuntimeError(f"Fig S5E expected 9 major x-ticks, found {len(maj_x)}: {maj_x}")
    xa = D.Axis("linear", [(maj_x[i], float(i * 30)) for i in range(9)])

    # Right major ticks: 1.0, 0.8, 0.6, 0.4, 0.2, 0.0
    r_ticks = find_inward_ticks(img, box, "right", depth_start=5, depth_end=28, min_depth=14, max_width=10)
    maj_y = sorted([p for p, d in r_ticks if rT + 50 < p < rB - 50])
    if len(maj_y) != 6:
        raise RuntimeError(f"Fig S5E expected 6 major y-ticks, found {len(maj_y)}: {maj_y}")
    ya = D.Axis("linear", [(maj_y[i], 1.0 - 0.2 * i) for i in range(6)])
    return xa, ya, box


def read_figS5e(img: np.ndarray, at: list[float], series: str) -> list[dict]:
    xa, ya, box = calibrate_figS5e(img)
    cL, rT, cR, rB = box
    r, g, b = img[:, :, 0].astype(np.float32), img[:, :, 1].astype(np.float32), img[:, :, 2].astype(np.float32)

    if series == "actin":
        mask = (g - r > 40) & (g - b > 20) & (g > 80)
    elif series == "bsa":
        mask = (r - b > 80) & (r > 170) & (g > 90)
    else:
        raise ValueError(series)

    # Mask legend box at top-right (x > 135 s, y > 0.78)
    mask[: int(round(ya.to_pixel(0.76))), int(round(xa.to_pixel(135.0))) :] = False

    # Special case: at pre-bleach t=-8 s and post-bleach t=0 s, orange BSA marker is drawn on top of green actin marker
    # Both share the exact same normalized pre-bleach (1.0) and t=0 (~0.0) reference points per SI Eq. (2).
    if series == "actin":
        orange = (r - b > 80) & (r > 170) & (g > 90)
        c_pre = int(round(xa.to_pixel(-2.0)))
        mask[:, :c_pre] |= orange[:, :c_pre]
        c_zero0 = int(round(xa.to_pixel(-2.0)))
        c_zero1 = int(round(xa.to_pixel(2.0)))
        r_zero0 = int(round(ya.to_pixel(0.05)))
        r_zero1 = int(round(ya.to_pixel(-0.05)))
        mask[r_zero0:r_zero1, c_zero0:c_zero1] |= orange[r_zero0:r_zero1, c_zero0:c_zero1]

    # Separate marker blobs from thin steep drop line at t=-8 s and t=0 s via morphological opening
    yy, xx = np.ogrid[-4:5, -4:5]
    disk_m = (xx * xx + yy * yy) <= 16
    opened = ndimage.binary_opening(mask, structure=disk_m)

    out = []
    for xv in at:
        # Use opened marker mask near the steep bleach transition (t <= 4 s) to avoid connecting line pull
        m_use = opened if xv <= 4.0 else mask
        res = _trace_mask_at(m_use, xa, ya, box, [xv], half_w=4 if xv <= 4.0 else 5, fallback_w=14, inset=4, extra_px=2.0)
        out.extend(res)
    return out


# ---------------------------------------------------------------------------
# Figure S6C (page 31, 600 dpi): Phalloidin fluorescence during filament winding
# ---------------------------------------------------------------------------

def calibrate_figS6c(img: np.ndarray):
    box = find_frame(img, dark_thr=140, sat_thr=40, min_h_frac=0.55, min_v_frac=0.55)
    cL, rT, cR, rB = box
    w, h = cR - cL, rB - rT
    # X ticks: 0 (cL), 600, 1200, 1800, 2400 (cR)
    b_ticks = find_inward_ticks(img, box, "bottom", depth_start=6, depth_end=26, min_depth=10, max_width=14)
    x_marks = [(cL, 0.0), (cR, 2400.0)]
    for val in (600.0, 1200.0, 1800.0):
        ep = cL + w * val / 2400.0
        near = [p for p, _ in b_ticks if abs(p - ep) <= 14]
        if near:
            x_marks.append((min(near, key=lambda p: abs(p - ep)), val))
    xa = D.Axis("linear", sorted(x_marks))

    # Y ticks on right spine: 12 (rT) down to 3 (rB), major every 1 unit
    r_ticks = find_inward_ticks(img, box, "right", depth_start=6, depth_end=30, min_depth=18, max_width=14)
    y_marks = [(rT, 12.0), (rB, 3.0)]
    for val in range(4, 12):
        ep = rB - h * (val - 3.0) / 9.0
        near = [p for p, _ in r_ticks if abs(p - ep) <= 14]
        if near:
            y_marks.append((min(near, key=lambda p: abs(p - ep)), float(val)))
    ya = D.Axis("linear", sorted(y_marks))
    return xa, ya, box


def read_figS6c(img: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, box = calibrate_figS6c(img)
    cL, rT, cR, rB = box
    gray, sat = _gray_sat(img)
    mask = (gray < 110) & (sat < 40)
    mask[: int(rT) + 25, :] = False
    mask[int(rB) - 25 :, :] = False
    mask[:, : int(cL) + 18] = False
    mask[:, int(cR) - 25 :] = False
    return _trace_mask_at(mask, xa, ya, box, at, half_w=3, fallback_w=10, inset=3, extra_px=2.0)
