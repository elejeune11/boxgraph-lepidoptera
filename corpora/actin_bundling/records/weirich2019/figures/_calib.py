"""Shared axis-calibration, curve-sampling, marker, and error-bar helpers for weirich2019."""

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
    """Fast connected-component extraction using scipy.ndimage."""
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


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 2.5) -> Axis:
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
        rs = runs_1d(sub[r], gap=4)
        if rs and max(b - a + 1 for a, b in rs) >= min_h:
            h_idx.append(int(r) + r0)

    v_idx = []
    for c in v_cand:
        rs = runs_1d(sub[:, c], gap=4)
        if rs and max(b - a + 1 for a, b in rs) >= min_v:
            v_idx.append(int(c) + c0)

    return dark, gray, clusters(np.array(h_idx, dtype=int), gap=3), clusters(np.array(v_idx, dtype=int), gap=3)


def find_ticks(
    mask: np.ndarray,
    spine: tuple[float, int, int],
    span: tuple[int, int],
    orientation: str,
    side: int,
    band: tuple[int, int] = (5, 20),
    frac: float = 0.55,
    max_width: int = 22,
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


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 5.0) -> list[float]:
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


def match_ticks_values(
    candidates: list[float],
    values: list[float],
    scale: str = "linear",
    tol: float = 6.0,
    increasing_pixel: bool = True,
) -> list[tuple[float, float]]:
    """Match a list of known tick values (sorted ascending) to pixel candidates.

    Uses every pair of candidates as anchor for (values[0], values[-1]) and picks
    the match minimizing max pixel residual.
    """
    c = sorted(candidates)
    us = np.array([math.log10(v) if scale == "log10" else float(v) for v in values], dtype=float)
    best_marks = None
    best_err = float("inf")
    n = len(values)
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            p0, p1 = (c[i], c[j]) if increasing_pixel else (c[j], c[i])
            if abs(p1 - p0) < 20.0:
                continue
            slope = (p1 - p0) / (us[-1] - us[0])
            pts = []
            ok = True
            for u in us:
                target = p0 + slope * (u - us[0])
                near = min(c, key=lambda x: abs(x - target))
                if abs(near - target) > tol:
                    ok = False
                    break
                pts.append(near)
            if ok and len(set(pts)) == n:
                ps = np.array(pts, dtype=float)
                b, a = np.polyfit(ps, us, 1)
                err = float(np.max(np.abs(ps - (us - a) / b)))
                if err < best_err:
                    best_err = err
                    best_marks = [(pts[k], float(values[k])) for k in range(n)]
    if best_marks is None:
        raise RuntimeError(f"Could not match {values} ({scale}) from candidates {c}")
    return best_marks


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
            if abs(nm - r_mid) <= max(45.0, 2.5 * thick):
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
    """Sample y(x) on an isolated dashed/dotted curve mask using weighted local polynomial regression."""
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
        seg = max(rs, key=lambda t: t[1] - t[0])
        col_xs.append(float(cc) + 0.5)
        col_ys.append(0.5 * (seg[0] + seg[1]) + 0.5 + row_min)
        col_thicks.append(float(seg[1] - seg[0] + 1))

    if len(col_xs) < 4:
        return {"x": float(x_val), "y": None, "found": False, "note": f"too few curve pixels near col {c}"}

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
    gap_unc = 0.0 if dist_nearest <= 2.0 else min(5.0, 0.05 * dist_nearest)
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
            f"{note_prefix}local poly fit (n={len(xs)}, d={dist_nearest:.0f}px, rms={fit_rms:.2f}px), "
            f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
        ),
    }


def make_point_reading(
    c: float,
    r: float,
    x_axis: Axis,
    y_axis: Axis,
    half_px_y: float = 3.0,
    note: str = "",
    error_bars: tuple[float, float, float, float] | None = None,
) -> dict:
    """Convert pixel center (c, r) to a calibrated reading dict."""
    x_val = float(x_axis.to_value(c))
    y_val = float(y_axis.to_value(r))
    h = float(half_px_y) + y_axis.residual_px
    y_a = float(y_axis.to_value(r - h))
    y_b = float(y_axis.to_value(r + h))
    out = {
        "x": x_val,
        "y": y_val,
        "lo": min(y_a, y_b),
        "hi": max(y_a, y_b),
        "px": [float(c), float(r)],
        "px_lo": float(r - h),
        "px_hi": float(r + h),
        "note": note,
    }
    if error_bars is not None:
        bar_lo, bar_hi, bar_lo_outer, bar_hi_outer = error_bars
        out["bar_lo"] = bar_lo
        out["bar_hi"] = bar_hi
        out["bar_lo_outer"] = bar_lo_outer
        out["bar_hi_outer"] = bar_hi_outer
    return out


def detect_template_mp(
    mask: np.ndarray,
    tpl: np.ndarray,
    min_full_frac: float = 0.60,
    min_res_frac: float = 0.40,
    min_dist: int = 16,
    refine_r: int = 5,
    max_hits: int = 200,
    stride: int = 2,
) -> list[tuple[float, float]]:
    """Fast matching pursuit for overlapping markers using 2x coarse FFT + full-res local refinement.

    Returns list of (row_center, col_center) in continuous coordinates.
    """
    from scipy import signal

    H, W = mask.shape
    th, tw = tpl.shape
    hr, hc = th // 2, tw // 2

    if stride > 1 and H * W > 600_000 and min(th, tw) >= 15:
        m_sub = mask[::stride, ::stride]
        t_sub = tpl[::stride, ::stride]
        coarse = detect_template_mp(
            m_sub,
            t_sub,
            min_full_frac=min_full_frac,
            min_res_frac=min_res_frac,
            min_dist=max(2, min_dist // stride),
            refine_r=max(2, refine_r // stride),
            max_hits=max_hits,
            stride=1,
        )
        # Refine every coarse center on the full-resolution mask inside a small +-refine_r patch
        tpl_f = tpl.astype(float)
        tpl_rev = tpl_f[::-1, ::-1]
        refined: list[tuple[float, float]] = []
        pad = refine_r + stride
        for rc, cc in coarse:
            r_est = int(round((rc - 0.5) * stride))
            c_est = int(round((cc - 0.5) * stride))
            pr0, pr1 = max(0, r_est - hr - pad), min(H, r_est + hr + pad + 1)
            pc0, pc1 = max(0, c_est - hc - pad), min(W, c_est + hc + pad + 1)
            patch = mask[pr0:pr1, pc0:pc1].astype(float)
            loc = signal.fftconvolve(patch, tpl_rev, mode="same")
            # Search only within +-pad of (r_est - pr0, c_est - pc0)
            rr0, rr1 = max(0, r_est - pr0 - pad), min(loc.shape[0], r_est - pr0 + pad + 1)
            cc0, cc1 = max(0, c_est - pc0 - pad), min(loc.shape[1], c_est - pc0 + pad + 1)
            sub = loc[rr0:rr1, cc0:cc1]
            dr, dc = np.unravel_index(np.argmax(sub), sub.shape)
            refined.append((float(pr0 + rr0 + int(dr)) + 0.5, float(pc0 + cc0 + int(dc)) + 0.5))
        return refined

    tpl_f = tpl.astype(float)
    tpl_area = float(tpl_f.sum())
    tpl_rev = tpl_f[::-1, ::-1]

    res = mask.astype(float).copy()
    full_score = signal.fftconvolve(res, tpl_rev, mode="same")
    score = full_score.copy()
    valid = full_score >= (min_full_frac * tpl_area)

    found: list[tuple[float, float]] = []
    for _ in range(max_hits):
        masked = score * valid
        idx = int(np.argmax(masked))
        r, c = divmod(idx, W)
        if masked[r, c] < min_res_frac * tpl_area:
            break
        r0, r1 = max(0, r - refine_r), min(H, r + refine_r + 1)
        c0, c1 = max(0, c - refine_r), min(W, c + refine_r + 1)
        sub = full_score[r0:r1, c0:c1]
        dr, dc = np.unravel_index(np.argmax(sub), sub.shape)
        rr, cc = r0 + int(dr), c0 + int(dc)
        found.append((float(rr) + 0.5, float(cc) + 0.5))

        valid[max(0, rr - min_dist) : min(H, rr + min_dist + 1), max(0, cc - min_dist) : min(W, cc + min_dist + 1)] = False

        pr0, pr1 = max(0, rr - hr), min(H, rr - hr + th)
        pc0, pc1 = max(0, cc - hc), min(W, cc - hc + tw)
        tr0, tr1 = pr0 - (rr - hr), pr1 - (rr - hr)
        tc0, tc1 = pc0 - (cc - hc), pc1 - (cc - hc)

        old_patch = res[pr0:pr1, pc0:pc1].copy()
        new_patch = old_patch * (1.0 - tpl_f[tr0:tr1, tc0:tc1])
        delta = old_patch - new_patch
        res[pr0:pr1, pc0:pc1] = new_patch

        conv = signal.fftconvolve(delta, tpl_rev, mode="full")
        sr0_raw, sc0_raw = pr0 - (th - 1 - hr), pc0 - (tw - 1 - hc)
        sr1_raw, sc1_raw = sr0_raw + conv.shape[0], sc0_raw + conv.shape[1]
        sr0, sr1 = max(0, sr0_raw), min(H, sr1_raw)
        sc0, sc1 = max(0, sc0_raw), min(W, sc1_raw)
        score[sr0:sr1, sc0:sc1] -= conv[sr0 - sr0_raw : sr1 - sr0_raw, sc0 - sc0_raw : sc1 - sc0_raw]

    return found

