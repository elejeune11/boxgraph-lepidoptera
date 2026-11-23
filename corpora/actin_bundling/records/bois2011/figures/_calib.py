"""Shared calibration and curve-reading helpers for Bois et al. (2011) figures."""

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


def find_frame(image: np.ndarray, row_frac: float = 0.5, col_frac: float = 0.45, thresh: float = 128.0):
    """Detect long horizontal and vertical plot-box frame lines."""
    gray = image.astype(float).mean(axis=2)
    dark = gray < thresh
    h_lines = clusters(np.where(dark.mean(axis=1) > row_frac)[0])
    v_lines = clusters(np.where(dark.mean(axis=0) > col_frac)[0])
    return dark, gray, h_lines, v_lines


def find_inward_ticks(
    dark: np.ndarray,
    h_top: tuple[float, int, int],
    h_bot: tuple[float, int, int],
    v_left: tuple[float, int, int],
    v_right: tuple[float, int, int],
    band: tuple[int, int] = (4, 13),
    margin: int = 15,
    max_width: int = 12,
    frac: float = 0.7,
):
    """Find inward-pointing ticks along all four borders of a plot box."""
    _, _, rt1 = h_top
    _, rb0, _ = h_bot
    _, _, cl1 = v_left
    _, cr0, _ = v_right
    b0, b1 = band

    def filt(cs):
        return [c[0] for c in cs if (c[2] - c[1] + 1) <= max_width]

    xt_top = filt(clusters(np.where(dark[rt1 + b0 : rt1 + b1, cl1 + margin : cr0 - margin].mean(axis=0) > frac)[0] + cl1 + margin))
    xt_bot = filt(clusters(np.where(dark[rb0 - b1 : rb0 - b0, cl1 + margin : cr0 - margin].mean(axis=0) > frac)[0] + cl1 + margin))
    yt_left = filt(clusters(np.where(dark[rt1 + margin : rb0 - margin, cl1 + b0 : cl1 + b1].mean(axis=1) > frac)[0] + rt1 + margin))
    yt_right = filt(clusters(np.where(dark[rt1 + margin : rb0 - margin, cr0 - b1 : cr0 - b0].mean(axis=1) > frac)[0] + rt1 + margin))
    return xt_top, xt_bot, yt_left, yt_right


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 2.0) -> list[float]:
    """Select the n_expected most evenly spaced tick positions from candidates."""
    c = sorted(candidates)
    if len(c) == n_expected:
        return c
    best = None
    best_err = float("inf")
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = (c[j] - c[i]) / (n_expected - 1)
            if step < 10.0:
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
                # linear fit max residual
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


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 1.5) -> Axis:
    """Build an Axis and assert its calibration residual is within tolerance."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_residual_px:
        raise RuntimeError(f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px: {marks}")
    return ax


def runs_1d(mask_1d: np.ndarray, gap: int = 1) -> list[tuple[int, int]]:
    """Return (start, end_inclusive) runs of True in a 1D boolean array."""
    idx = np.nonzero(mask_1d)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap + 1)[0] + 1
    return [(int(g[0]), int(g[-1])) for g in np.split(idx, splits)]


def sample_solid_curve_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    row_min: int,
    row_max: int,
    pick_fn=None,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on a solid curve mask at x_val within [row_min, row_max].

    pick_fn(runs, col) selects the target (r0, r1) run if multiple runs exist in that column;
    by default requires a single run or picks the thickest run if none given.
    """
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c = max(0, min(W - 1, c))
    sub = mask[row_min : row_max + 1, c]
    rs = [(r0 + row_min, r1 + row_min) for r0, r1 in runs_1d(sub)]
    if not rs:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no curve pixel at col {c}"}
    if pick_fn is not None:
        seg = pick_fn(rs, c)
    else:
        seg = rs[0] if len(rs) == 1 else max(rs, key=lambda t: t[1] - t[0])
    if seg is None:
        return {"x": float(x_val), "y": None, "found": False, "note": f"no matching run at col {c}"}

    r0, r1 = seg
    r_mid = 0.5 * (r0 + r1) + 0.5
    thick = float(r1 - r0 + 1)

    # Estimate local slope dr/dc using columns c-3 and c+3
    slope = 0.0
    pts = []
    for dc in (-4, -2, 0, 2, 4):
        cc = c + dc
        if 0 <= cc < W:
            s_rs = [(a + row_min, b + row_min) for a, b in runs_1d(mask[row_min : row_max + 1, cc])]
            if s_rs:
                near = min(s_rs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - r_mid))
                if abs(0.5 * (near[0] + near[1]) + 0.5 - r_mid) <= max(25.0, 2.0 * thick):
                    pts.append((float(cc) + 0.5, 0.5 * (near[0] + near[1]) + 0.5))
    if len(pts) >= 2:
        xs_p = np.array([p[0] for p in pts])
        ys_p = np.array([p[1] for p in pts])
        slope = float(np.polyfit(xs_p, ys_p, 1)[0])

    # Interpolate subpixel offset from integer column center (c + 0.5) to exact c_f
    r_est = r_mid + slope * (c_f - (c + 0.5))
    # Honest 95% half-width in pixels: half stroke vertical height + quantization + calibration misfit
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
    """Sample y(x) on an isolated dashed/dotted curve mask using a local polynomial fit
    to all curve pixels within +/- max_search_dx columns, eliminating bias from slanted dash tips."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c_lo = max(0, c - max_search_dx)
    c_hi = min(W - 1, c + max_search_dx)
    sub = mask[row_min : row_max + 1, c_lo : c_hi + 1]
    rr, cc = np.nonzero(sub)
    if len(rr) < 15:
        return {"x": float(x_val), "y": None, "found": False, "note": f"too few dash/dot pixels near col {c}"}

    xs = cc.astype(float) + c_lo + 0.5
    ys = rr.astype(float) + row_min + 0.5

    # Ensure pixels exist on both sides of c_f (or very close to c_f) so we don't extrapolate wildly
    has_left = bool(np.any(xs <= c_f + 2.0))
    has_right = bool(np.any(xs >= c_f - 2.0))
    if not (has_left and has_right):
        return {"x": float(x_val), "y": None, "found": False, "note": f"one-sided dashes near col {c}"}

    # Gaussian weights centered at c_f with sigma = 0.5 * max_search_dx
    sigma = max(18.0, 0.55 * max_search_dx)
    w = np.exp(-0.5 * ((xs - c_f) / sigma) ** 2)

    # Reject any outlier pixels far from a robust linear fit (in case another stray speck entered the window)
    p1 = np.poly1d(np.polyfit(xs, ys, 1, w=w))
    resid1 = np.abs(ys - p1(xs))
    keep = resid1 <= max(35.0, 3.0 * float(np.median(resid1)) + 15.0)
    xs, ys, w = xs[keep], ys[keep], w[keep]

    # Fit local quadratic (and linear for curvature check)
    p_lin = np.poly1d(np.polyfit(xs, ys, 1, w=w))
    p_quad = np.poly1d(np.polyfit(xs, ys, 2, w=w))
    r_est = float(p_quad(c_f))
    slope = float(np.polyder(p_quad)(c_f))

    # Estimate vertical stroke thickness from 95% spread of residuals around p_quad
    resid = np.abs(ys - p_quad(xs))
    half_thick = max(float(np.percentile(resid, 95)), 4.0)
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
        "note": f"{note_prefix}n_px={len(xs)}, h_thick={half_thick:.1f}px, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
    }


def sample_horizontal_x_of_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    y_val: float,
    col_min: int,
    col_max: int,
    pick_fn=None,
    note_prefix: str = "",
) -> dict:
    """Sample horizontal coordinate x(y) at vertical value y_val (returns datum x=y_val, y=x_val)."""
    r_f = float(y_axis.to_pixel(y_val))
    r = int(round(r_f))
    H, W = mask.shape
    r = max(0, min(H - 1, r))
    cs = [(c0 + col_min, c1 + col_min) for c0, c1 in runs_1d(mask[r, col_min : col_max + 1])]
    if not cs:
        return {"x": float(y_val), "y": None, "found": False, "note": f"no curve pixel at row {r}"}
    if pick_fn is not None:
        seg = pick_fn(cs, r)
    else:
        seg = cs[0] if len(cs) == 1 else max(cs, key=lambda t: t[1] - t[0])
    if seg is None:
        return {"x": float(y_val), "y": None, "found": False, "note": f"no matching run at row {r}"}

    c0, c1 = seg
    c_mid = 0.5 * (c0 + c1) + 0.5
    thick = float(c1 - c0 + 1)

    slope = 0.0
    pts = []
    for dr in (-4, -2, 0, 2, 4):
        rr = r + dr
        if 0 <= rr < H:
            s_cs = [(a + col_min, b + col_min) for a, b in runs_1d(mask[rr, col_min : col_max + 1])]
            if s_cs:
                near = min(s_cs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - c_mid))
                if abs(0.5 * (near[0] + near[1]) + 0.5 - c_mid) <= max(25.0, 2.0 * thick):
                    pts.append((float(rr) + 0.5, 0.5 * (near[0] + near[1]) + 0.5))
    if len(pts) >= 2:
        rs_p = np.array([p[0] for p in pts])
        cs_p = np.array([p[1] for p in pts])
        slope = float(np.polyfit(rs_p, cs_p, 1)[0])

    c_est = c_mid + slope * (r_f - (r + 0.5))
    half_px = max(0.5 * thick, 1.5) + 0.5 + x_axis.residual_px + abs(slope) * (0.5 + y_axis.residual_px)
    val = float(x_axis.to_value(c_est))
    v_a = float(x_axis.to_value(c_est - half_px))
    v_b = float(x_axis.to_value(c_est + half_px))
    lo, hi = min(v_a, v_b), max(v_a, v_b)
    return {
        "x": float(y_val),
        "y": val,
        "lo": lo,
        "hi": hi,
        "px": [float(c_est), r_f],
        "note": f"{note_prefix}cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, thick={thick:.0f}px",
    }
