"""Shared axis-calibration, curve-sampling, marker, and error-bar helpers for weirich2017."""

from __future__ import annotations

import math
import numpy as np
from boxgraph.digitize import Axis


def clusters(indices: np.ndarray, gap: int = 3) -> list[tuple[float, int, int]]:
    """Group sorted integer indices into contiguous clusters (center, start, end)."""
    if len(indices) == 0:
        return []
    splits = np.where(np.diff(indices) > gap)[0] + 1
    out = []
    for g in np.split(indices, splits):
        out.append((float(g.mean()) + 0.5, int(g[0]), int(g[-1])))
    return out


def runs_1d(mask_1d: np.ndarray, gap: int = 1) -> list[tuple[int, int]]:
    """Return (start, end_inclusive) runs of True in a 1D boolean array."""
    idx = np.nonzero(mask_1d)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap + 1)[0] + 1
    return [(int(g[0]), int(g[-1])) for g in np.split(idx, splits)]


def find_blobs(mask: np.ndarray, min_area: int = 10, max_area: int | None = None) -> list[dict]:
    """Fast connected-component extraction using ndimage.find_objects."""
    from scipy import ndimage

    lab, _ = ndimage.label(mask)
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        if sl is None:
            continue
        sub = lab[sl] == i
        area = int(sub.sum())
        if area < min_area or (max_area is not None and area > max_area):
            continue
        rr, cc = np.nonzero(sub)
        r0, c0 = sl[0].start, sl[1].start
        r1, c1 = sl[0].stop - 1, sl[1].stop - 1
        out.append(
            {
                "label": i,
                "area": area,
                "r": float(rr.mean()) + r0 + 0.5,
                "c": float(cc.mean()) + c0 + 0.5,
                "r0": r0,
                "r1": r1,
                "c0": c0,
                "c1": c1,
                "h": r1 - r0 + 1,
                "w": c1 - c0 + 1,
            }
        )
    return out


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 2.0) -> Axis:
    """Build an Axis and verify its calibration residual is within tolerance."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_residual_px:
        raise RuntimeError(
            f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px: {marks}"
        )
    return ax


def find_box_spines(
    image: np.ndarray,
    thresh: float = 128.0,
    row_frac: float = 0.45,
    col_frac: float = 0.45,
    row_range: tuple[int, int] | None = None,
    col_range: tuple[int, int] | None = None,
    dark_mask: np.ndarray | None = None,
):
    """Find long horizontal and vertical spines in an image or subregion.

    Returns (dark, gray, h_lines, v_lines) where each line is (center, start, end).
    """
    gray = image.astype(float).mean(axis=2)
    dark = dark_mask if dark_mask is not None else (gray < thresh)
    H, W = dark.shape
    r0, r1 = row_range if row_range is not None else (0, H)
    c0, c1 = col_range if col_range is not None else (0, W)
    sub = dark[r0:r1, c0:c1]

    min_h = row_frac * sub.shape[1]
    min_v = col_frac * sub.shape[0]
    h_cand = np.where(sub.sum(axis=1) >= min_h)[0]
    v_cand = np.where(sub.sum(axis=0) >= min_v)[0]

    h_idx = []
    for r in h_cand:
        rs = runs_1d(sub[r], gap=3)
        if rs and max(b - a + 1 for a, b in rs) >= min_h:
            h_idx.append(int(r) + r0)

    v_idx = []
    for c in v_cand:
        rs = runs_1d(sub[:, c], gap=3)
        if rs and max(b - a + 1 for a, b in rs) >= min_v:
            v_idx.append(int(c) + c0)

    return dark, gray, clusters(np.array(h_idx, dtype=int), gap=3), clusters(np.array(v_idx, dtype=int), gap=3)


def find_ticks(
    mask: np.ndarray,
    spine: tuple[float, int, int],
    span: tuple[int, int],
    orientation: str,
    side: int,
    band: tuple[int, int] = (6, 22),
    frac: float = 0.55,
    max_width: int = 18,
    gap: int = 2,
) -> list[float]:
    """Find tick positions protruding from a spine in direction `side` (-1 or +1).

    orientation='x': ticks stick vertically off horizontal `spine`; `span` = (c_min, c_max).
                     side=-1 looks above (inward from bottom spine), side=+1 looks below.
    orientation='y': ticks stick horizontally off vertical `spine`; `span` = (r_min, r_max).
                     side=+1 looks right (inward from left spine), side=-1 looks left.
    """
    _, s0, s1 = spine
    a, b = span
    b0, b1 = band
    H, W = mask.shape
    if orientation == "y":
        if side < 0:
            c0, c1 = max(0, s0 - b1), max(0, s0 - b0)
        else:
            c0, c1 = min(W, s1 + 1 + b0), min(W, s1 + 1 + b1)
        prof = mask[max(0, a) : min(H, b), c0:c1].mean(axis=1)
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=gap)
    elif orientation == "x":
        if side > 0:
            r0, r1 = min(H, s1 + 1 + b0), min(H, s1 + 1 + b1)
        else:
            r0, r1 = max(0, s0 - b1), max(0, s0 - b0)
        prof = mask[r0:r1, max(0, a) : min(W, b)].mean(axis=0)
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=gap)
    else:
        raise ValueError(f"Unknown orientation {orientation!r}")
    return [c[0] for c in cs if (c[2] - c[1] + 1) <= max_width]


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 4.0) -> list[float]:
    """Select the n_expected most evenly spaced tick positions from candidates."""
    c = sorted(candidates)
    if len(c) == n_expected:
        return c
    best = None
    best_err = float("inf")
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = (c[j] - c[i]) / (n_expected - 1)
            if step < 15.0:
                continue
            pts = [c[i]]
            ok = True
            for k in range(1, n_expected):
                target = c[i] + k * step
                near = min(c, key=lambda x: abs(x - target))
                if abs(near - target) > tol:
                    ok = False
                    break
                pts.append(near)
            if ok and len(set(pts)) == n_expected:
                xs = np.arange(n_expected, dtype=float)
                ps = np.array(pts, dtype=float)
                m, b = np.polyfit(xs, ps, 1)
                err = float(np.max(np.abs(ps - (m * xs + b))))
                if err < best_err:
                    best_err = err
                    best = pts
    if best is None:
        raise RuntimeError(f"Could not match {n_expected} regular ticks from {c}")
    return best


def sample_solid_curve_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    row_min: int,
    row_max: int,
    pick_fn=None,
    search_dx: int = 0,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on a solid curve mask at x_val within [row_min, row_max]."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c = max(0, min(W - 1, c))

    def get_runs(col):
        if col < 0 or col >= W:
            return []
        return [(r0 + row_min, r1 + row_min) for r0, r1 in runs_1d(mask[row_min : row_max + 1, col])]

    rs = get_runs(c)
    used_c = c
    if not rs and search_dx > 0:
        for d in range(1, search_dx + 1):
            for cc in (c - d, c + d):
                rs = get_runs(cc)
                if rs:
                    used_c = cc
                    break
            if rs:
                break
    if not rs:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no curve pixel at col {c}"}

    if pick_fn is not None:
        seg = pick_fn(rs, used_c)
    else:
        seg = rs[0] if len(rs) == 1 else max(rs, key=lambda t: t[1] - t[0])
    if seg is None:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no matching run at col {used_c}"}

    r0, r1 = seg
    r_mid = 0.5 * (r0 + r1) + 0.5
    thick = float(r1 - r0 + 1)

    pts = []
    for dc in (-5, -3, -1, 0, 1, 3, 5):
        cc = used_c + dc
        s_rs = get_runs(cc)
        if s_rs:
            near = min(s_rs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - r_mid))
            nm = 0.5 * (near[0] + near[1]) + 0.5
            if abs(nm - r_mid) <= max(40.0, 2.5 * thick):
                pts.append((float(cc) + 0.5, nm))
    slope = 0.0
    if len(pts) >= 2:
        xs_p = np.array([p[0] for p in pts])
        ys_p = np.array([p[1] for p in pts])
        if np.ptp(xs_p) > 0:
            slope = float(np.polyfit(xs_p, ys_p, 1)[0])

    r_est = r_mid + slope * (c_f - (used_c + 0.5))
    half_px = max(0.5 * thick, 1.5) + 0.5 + y_axis.residual_px + abs(slope) * (0.5 + x_axis.residual_px)
    y_val = float(y_axis.to_value(r_est))
    y_a = float(y_axis.to_value(r_est - half_px))
    y_b = float(y_axis.to_value(r_est + half_px))
    lo, hi = min(y_a, y_b), max(y_a, y_b)
    return {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [c_f, float(r_est)],
        "px_lo": float(r_est - half_px),
        "px_hi": float(r_est + half_px),
        "note": f"{note_prefix}cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, thick={thick:.0f}px",
    }


def sample_dashed_curve_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    row_min: int,
    row_max: int,
    max_search_dx: int = 120,
    sigma_dx: float = 45.0,
    deg: int = 2,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on an isolated dashed curve mask using weighted local polynomial regression."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c_lo = max(0, c - max_search_dx)
    c_hi = min(W - 1, c + max_search_dx)

    col_xs = []
    col_ys = []
    col_thicks = []
    for cc in range(c_lo, c_hi + 1):
        rs = runs_1d(mask[row_min : row_max + 1, cc], gap=2)
        if not rs:
            continue
        # if multiple runs in column, take the thickest or median-closest (mask should be isolated)
        seg = max(rs, key=lambda t: t[1] - t[0])
        col_xs.append(float(cc) + 0.5)
        col_ys.append(0.5 * (seg[0] + seg[1]) + 0.5 + row_min)
        col_thicks.append(float(seg[1] - seg[0] + 1))

    if len(col_xs) < 4:
        return {"x": float(x_val), "y": None, "found": False, "note": f"too few dashed pixels near col {c}"}

    xs = np.array(col_xs) - c_f
    ys = np.array(col_ys)
    w = np.exp(-0.5 * (xs / sigma_dx) ** 2)
    eff_deg = min(deg, 1 if (np.ptp(xs) < 25 or len(xs) < 10) else deg)
    coeffs = np.polyfit(xs, ys, eff_deg, w=np.sqrt(w))
    r_est = float(np.polyval(coeffs, 0.0))
    slope = float(np.polyval(np.polyder(coeffs), 0.0))
    resid = np.abs(ys - np.polyval(coeffs, xs))
    fit_rms = float(np.sqrt(np.average(resid**2, weights=w)))
    thick = float(np.median(col_thicks))

    dist_nearest = float(np.min(np.abs(xs)))
    gap_unc = 0.0 if dist_nearest <= 2.0 else min(4.0, 0.04 * dist_nearest)
    half_px = max(0.5 * thick, 1.5) + 0.5 + fit_rms + gap_unc + y_axis.residual_px + abs(slope) * (0.5 + x_axis.residual_px)
    y_val = float(y_axis.to_value(r_est))
    y_a = float(y_axis.to_value(r_est - half_px))
    y_b = float(y_axis.to_value(r_est + half_px))
    lo, hi = min(y_a, y_b), max(y_a, y_b)
    return {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [c_f, float(r_est)],
        "px_lo": float(r_est - half_px),
        "px_hi": float(r_est + half_px),
        "note": (
            f"{note_prefix}dashed fit (n={len(xs)}, d={dist_nearest:.0f}px, rms={fit_rms:.2f}px), "
            f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
        ),
    }


def trace_dashed_error_bar(
    stem_mask: np.ndarray,
    col_c: float,
    r_top_start: int,
    r_bot_start: int,
    r_min_bound: int,
    r_max_bound: int,
    y_axis: Axis,
    stem_half_w: int = 5,
    gap: int = 6,
) -> tuple[float, float, float, float, float, float]:
    """Trace vertical error bar stem upward from r_top_start and downward from r_bot_start.

    Returns (bar_lo, bar_hi, bar_lo_outer, bar_hi_outer, row_hi, row_lo) in y-value and pixel coords.
    """
    H, W = stem_mask.shape
    c0 = max(0, int(round(col_c)) - stem_half_w)
    c1 = min(W, int(round(col_c)) + stem_half_w + 1)

    # Walk upward from r_top_start
    r_up = r_top_start
    missed = 0
    for r in range(r_top_start, max(r_min_bound - 1, -1), -1):
        if stem_mask[r, c0:c1].any():
            r_up = r
            missed = 0
        else:
            missed += 1
            if missed > gap:
                break

    # Walk downward from r_bot_start
    r_dn = r_bot_start
    missed = 0
    for r in range(r_bot_start, min(r_max_bound + 1, H)):
        if stem_mask[r, c0:c1].any():
            r_dn = r
            missed = 0
        else:
            missed += 1
            if missed > gap:
                break

    # Refine cap center rows if horizontal cap exists near r_up / r_dn
    def refine_cap(r_end: int, direction: int) -> float:
        # Look within 8 px inward from r_end for widest horizontal run
        best_r = float(r_end) + 0.5
        best_w = 0
        c_wide0 = max(0, int(round(col_c)) - 22)
        c_wide1 = min(W, int(round(col_c)) + 23)
        cap_rows = []
        for dr in range(0, 10):
            rr = r_end + direction * dr
            if 0 <= rr < H:
                w = int(stem_mask[rr, c_wide0:c_wide1].sum())
                if w > best_w:
                    best_w = w
        if best_w >= 12:
            for dr in range(0, 10):
                rr = r_end + direction * dr
                if 0 <= rr < H and int(stem_mask[rr, c_wide0:c_wide1].sum()) >= 0.7 * best_w:
                    cap_rows.append(rr)
            if cap_rows:
                best_r = float(np.mean(cap_rows)) + 0.5
        return best_r

    row_hi = refine_cap(r_up, +1)
    row_lo = refine_cap(r_dn, -1)

    end_unc = 2.0 + y_axis.residual_px
    y_hi = float(y_axis.to_value(row_hi))
    y_lo = float(y_axis.to_value(row_lo))
    y_hi_out = float(y_axis.to_value(row_hi - end_unc))
    y_lo_out = float(y_axis.to_value(row_lo + end_unc))
    bar_lo, bar_hi = min(y_lo, y_hi), max(y_lo, y_hi)
    bar_lo_outer, bar_hi_outer = min(y_lo_out, y_hi_out), max(y_lo_out, y_hi_out)
    return bar_lo, bar_hi, bar_lo_outer, bar_hi_outer, row_hi, row_lo
