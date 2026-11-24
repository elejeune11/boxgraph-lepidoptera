"""Shared calibration and curve-reading helpers for Frohoff-Hulsmann et al. (2021) figures."""

from __future__ import annotations

import math
import numpy as np
from boxgraph.digitize import Axis


def clusters(indices: np.ndarray, gap: int = 3) -> list[tuple[float, int, int]]:
    """Group sorted integer indices into contiguous clusters (center, start, end_inclusive)."""
    if len(indices) == 0:
        return []
    splits = np.where(np.diff(indices) > gap)[0] + 1
    out = []
    for g in np.split(indices, splits):
        out.append((float(g.mean()) + 0.5, int(g[0]), int(g[-1])))
    return out


def find_main_box(
    image: np.ndarray,
    row_frac: float = 0.45,
    col_frac: float = 0.45,
    thresh: float = 128.0,
):
    """Detect the outermost plot frame (h_top, h_bot, v_left, v_right).

    Returns (dark, gray, h_top, h_bot, v_left, v_right) where each frame line is
    a tuple (center_px, start_idx, end_idx).
    """
    gray = image.astype(float).mean(axis=2)
    # Frame lines are neutral black/dark gray (low chroma)
    chroma = image.max(axis=2).astype(float) - image.min(axis=2).astype(float)
    dark = (gray < thresh) & (chroma < 40.0)
    h_lines = clusters(np.where(dark.mean(axis=1) > row_frac)[0])
    v_lines = clusters(np.where(dark.mean(axis=0) > col_frac)[0])
    if len(h_lines) < 2 or len(v_lines) < 2:
        raise RuntimeError(
            f"Could not detect rectangular plot box: h_lines={h_lines}, v_lines={v_lines}"
        )
    return dark, gray, h_lines[0], h_lines[-1], v_lines[0], v_lines[-1]


def find_outward_ticks(
    dark: np.ndarray,
    h_top: tuple[float, int, int],
    h_bot: tuple[float, int, int],
    v_left: tuple[float, int, int],
    v_right: tuple[float, int, int],
    band: tuple[int, int] = (3, 10),
    pad: int = 6,
    max_width: int = 12,
    frac: float = 0.65,
):
    """Detect outward-pointing matplotlib ticks along the bottom and left borders."""
    H, W = dark.shape
    _, rt0, _ = h_top
    _, _, rb1 = h_bot
    _, cl0, _ = v_left
    _, _, cr1 = v_right
    b0, b1 = band

    def filt(cs):
        return [c[0] for c in cs if (c[2] - c[1] + 1) <= max_width]

    c_start = max(0, cl0 - pad)
    c_end = min(W, cr1 + pad + 1)
    r_b0 = min(H, rb1 + b0)
    r_b1 = min(H, rb1 + b1)
    xt_bot = filt(
        clusters(np.where(dark[r_b0:r_b1, c_start:c_end].mean(axis=0) > frac)[0] + c_start)
    )

    r_start = max(0, rt0 - pad)
    r_end = min(H, rb1 + pad + 1)
    c_l0 = max(0, cl0 - b1)
    c_l1 = max(0, cl0 - b0)
    yt_left = filt(
        clusters(np.where(dark[r_start:r_end, c_l0:c_l1].mean(axis=1) > frac)[0] + r_start)
    )
    return xt_bot, yt_left


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 2.5) -> list[float]:
    """Select the n_expected most evenly spaced tick positions from candidates."""
    c = sorted(candidates)
    if len(c) == n_expected:
        return c
    if len(c) < n_expected:
        raise RuntimeError(f"Expected {n_expected} ticks, only found {len(c)}: {c}")
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


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 1.2) -> Axis:
    """Build an Axis and verify its calibration residual is small."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_residual_px:
        raise RuntimeError(
            f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px: {marks}"
        )
    return ax


def rgb_mask(
    image: np.ndarray,
    target_rgb: tuple[float, float, float],
    tol: float = 55.0,
    min_chroma: float = 25.0,
) -> np.ndarray:
    """Boolean mask of pixels close to target_rgb with minimum chromaticity."""
    arr = image.astype(float)
    t = np.array(target_rgb, dtype=float)
    d = np.sqrt(np.sum((arr - t) ** 2, axis=2))
    mask = d <= tol
    if min_chroma > 0.0:
        chroma = arr.max(axis=2) - arr.min(axis=2)
        mask &= chroma >= min_chroma
    return mask


def runs_1d(mask_1d: np.ndarray, gap: int = 1) -> list[tuple[int, int]]:
    """Return (start, end_inclusive) runs of True in a 1D boolean array."""
    idx = np.nonzero(mask_1d)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap + 1)[0] + 1
    return [(int(g[0]), int(g[-1])) for g in np.split(idx, splits)]


def mask_interior(
    mask: np.ndarray,
    h_top: tuple[float, int, int],
    h_bot: tuple[float, int, int],
    v_left: tuple[float, int, int],
    v_right: tuple[float, int, int],
    inset: int = 4,
) -> np.ndarray:
    """Zero out everything outside the plot frame interior."""
    out = np.zeros_like(mask, dtype=bool)
    r0 = h_top[2] + inset
    r1 = h_bot[1] - inset
    c0 = v_left[2] + inset
    c1 = v_right[1] - inset
    out[r0 : r1 + 1, c0 : c1 + 1] = mask[r0 : r1 + 1, c0 : c1 + 1]
    return out


def sample_solid_curve_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    row_min: int,
    row_max: int,
    pick_fn=None,
    search_dc: int = 2,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on a solid curve mask at x_val within [row_min, row_max]."""
    c_f = float(x_axis.to_pixel(x_val))
    c0 = int(round(c_f))
    H, W = mask.shape
    c0 = max(0, min(W - 1, c0))

    rs = []
    c_used = c0
    for dc in sorted(range(-search_dc, search_dc + 1), key=abs):
        cc = c0 + dc
        if 0 <= cc < W:
            sub = mask[row_min : row_max + 1, cc]
            cand = [(a + row_min, b + row_min) for a, b in runs_1d(sub)]
            if cand:
                rs = cand
                c_used = cc
                break
    if not rs:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no curve pixel near col {c0}"}

    if pick_fn is not None:
        seg = pick_fn(rs, c_used)
    else:
        seg = rs[0] if len(rs) == 1 else max(rs, key=lambda t: t[1] - t[0])
    if seg is None:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no matching run at col {c_used}"}

    r0, r1 = seg
    r_mid = 0.5 * (r0 + r1) + 0.5
    thick = float(r1 - r0 + 1)

    slope = 0.0
    pts = []
    for dc in (-4, -2, 0, 2, 4):
        cc = c_used + dc
        if 0 <= cc < W:
            s_rs = [(a + row_min, b + row_min) for a, b in runs_1d(mask[row_min : row_max + 1, cc])]
            if s_rs:
                near = min(s_rs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - r_mid))
                if abs(0.5 * (near[0] + near[1]) + 0.5 - r_mid) <= max(30.0, 2.5 * thick):
                    pts.append((float(cc) + 0.5, 0.5 * (near[0] + near[1]) + 0.5))
    if len(pts) >= 2:
        xs_p = np.array([p[0] for p in pts])
        ys_p = np.array([p[1] for p in pts])
        slope = float(np.polyfit(xs_p, ys_p, 1)[0])

    r_est = r_mid + slope * (c_f - (c_used + 0.5))
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
    max_search_dx: int = 55,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on an isolated dashed/dotted/dash-dot curve via local quadratic regression."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c_lo = max(0, c - max_search_dx)
    c_hi = min(W - 1, c + max_search_dx)
    sub = mask[row_min : row_max + 1, c_lo : c_hi + 1]
    rr, cc = np.nonzero(sub)
    if len(rr) < 12:
        return {"x": float(x_val), "y": None, "found": False, "note": f"too few dash pixels near col {c}"}

    xs = cc.astype(float) + c_lo + 0.5
    ys = rr.astype(float) + row_min + 0.5

    has_left = bool(np.any(xs <= c_f + 3.0))
    has_right = bool(np.any(xs >= c_f - 3.0))
    if not (has_left and has_right):
        return {"x": float(x_val), "y": None, "found": False, "note": f"one-sided dashes near col {c}"}

    sigma = max(18.0, 0.55 * max_search_dx)
    w = np.exp(-0.5 * ((xs - c_f) / sigma) ** 2)

    p1 = np.poly1d(np.polyfit(xs, ys, 1, w=w))
    resid1 = np.abs(ys - p1(xs))
    keep = resid1 <= max(30.0, 3.0 * float(np.median(resid1)) + 12.0)
    xs, ys, w = xs[keep], ys[keep], w[keep]
    if len(xs) < 10:
        return {"x": float(x_val), "y": None, "found": False, "note": f"too few inlier pixels near col {c}"}

    p_lin = np.poly1d(np.polyfit(xs, ys, 1, w=w))
    p_quad = np.poly1d(np.polyfit(xs, ys, 2, w=w))
    r_est = float(p_quad(c_f))
    slope = float(np.polyder(p_quad)(c_f))

    resid = np.abs(ys - p_quad(xs))
    half_thick = max(float(np.percentile(resid, 95)), 3.0)
    curv_diff = 0.25 * abs(r_est - float(p_lin(c_f)))

    half_px = (
        half_thick
        + 1.0
        + curv_diff
        + y_axis.residual_px
        + abs(slope) * (0.5 + x_axis.residual_px)
    )
    y_val = float(y_axis.to_value(r_est))
    y_a = float(y_axis.to_value(r_est - half_px))
    y_b = float(y_axis.to_value(r_est + half_px))
    lo, hi = min(y_a, y_b), max(y_a, y_b)
    return {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [c_f, r_est],
        "px_lo": float(r_est - half_px),
        "px_hi": float(r_est + half_px),
        "note": f"{note_prefix}cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, half={half_px:.1f}px",
    }


def sample_horizontal_x_of_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    y_val: float,
    col_min: int,
    col_max: int,
    pick_fn=None,
    search_dr: int = 2,
    note_prefix: str = "",
) -> dict:
    """Sample horizontal x-coordinate as a function of vertical y-coordinate.

    Returns a reading dict where `x` is the queried vertical-axis coordinate (`y_val`)
    and `y` is the measured horizontal-axis value, matching BoxGraph horizontal reading convention.
    """
    r_f = float(y_axis.to_pixel(y_val))
    r0 = int(round(r_f))
    H, W = mask.shape
    r0 = max(0, min(H - 1, r0))

    cs = []
    r_used = r0
    for dr in sorted(range(-search_dr, search_dr + 1), key=abs):
        rr = r0 + dr
        if 0 <= rr < H:
            sub = mask[rr, col_min : col_max + 1]
            cand = [(a + col_min, b + col_min) for a, b in runs_1d(sub)]
            if cand:
                cs = cand
                r_used = rr
                break
    if not cs:
        return {"x": float(y_val), "y": None, "found": False, "note": f"no curve pixel near row {r0}"}

    if pick_fn is not None:
        seg = pick_fn(cs, r_used)
    else:
        seg = cs[0] if len(cs) == 1 else max(cs, key=lambda t: t[1] - t[0])
    if seg is None:
        return {"x": float(y_val), "y": None, "found": False, "note": f"no matching run at row {r_used}"}

    c0, c1 = seg
    c_mid = 0.5 * (c0 + c1) + 0.5
    thick = float(c1 - c0 + 1)

    slope = 0.0
    pts = []
    for dr in (-4, -2, 0, 2, 4):
        rr = r_used + dr
        if 0 <= rr < H:
            s_cs = [(a + col_min, b + col_min) for a, b in runs_1d(mask[rr, col_min : col_max + 1])]
            if s_cs:
                near = min(s_cs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - c_mid))
                if abs(0.5 * (near[0] + near[1]) + 0.5 - c_mid) <= max(30.0, 2.5 * thick):
                    pts.append((float(rr) + 0.5, 0.5 * (near[0] + near[1]) + 0.5))
    if len(pts) >= 2:
        rs_p = np.array([p[0] for p in pts])
        cs_p = np.array([p[1] for p in pts])
        slope = float(np.polyfit(rs_p, cs_p, 1)[0])

    c_est = c_mid + slope * (r_f - (r_used + 0.5))
    half_px = max(0.5 * thick, 1.5) + 0.5 + x_axis.residual_px + abs(slope) * (0.5 + y_axis.residual_px)
    x_measured = float(x_axis.to_value(c_est))
    x_a = float(x_axis.to_value(c_est - half_px))
    x_b = float(x_axis.to_value(c_est + half_px))
    lo, hi = min(x_a, x_b), max(x_a, x_b)
    return {
        "x": float(y_val),
        "y": x_measured,
        "lo": lo,
        "hi": hi,
        "px": [float(c_est), r_f],
        "note": f"{note_prefix}cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, thick={thick:.0f}px",
    }
