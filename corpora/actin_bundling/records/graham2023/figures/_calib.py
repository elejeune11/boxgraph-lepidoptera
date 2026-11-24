"""Shared calibration, curve, bar, boxplot, and heatmap reading helpers for graham2023."""

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


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 2.0) -> Axis:
    """Build an Axis and assert its calibration residual is within tolerance."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_residual_px:
        raise RuntimeError(
            f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px: {marks}"
        )
    return ax


def find_l_spines(
    image: np.ndarray,
    thresh: float = 128.0,
    row_frac: float = 0.35,
    col_frac: float = 0.35,
    row_range: tuple[int, int] | None = None,
    col_range: tuple[int, int] | None = None,
):
    """Find long dark horizontal and vertical axis lines (spines) in an image or subregion.

    Returns (dark, gray, h_lines, v_lines) where each line is (center, start, end) in image coords.
    """
    gray = image.astype(float).mean(axis=2)
    dark = gray < thresh
    H, W = dark.shape
    r0, r1 = row_range if row_range is not None else (0, H)
    c0, c1 = col_range if col_range is not None else (0, W)
    sub = dark[r0:r1, c0:c1]

    # Fast vectorized pre-filter: a contiguous run of length >= frac*N requires total sum >= frac*N
    min_h = row_frac * sub.shape[1]
    min_v = col_frac * sub.shape[0]
    h_cand = np.where(sub.sum(axis=1) >= min_h)[0]
    v_cand = np.where(sub.sum(axis=0) >= min_v)[0]
    h_idx = []
    for r in h_cand:
        rs = runs_1d(sub[r], gap=2)
        if rs and max(b - a + 1 for a, b in rs) >= min_h:
            h_idx.append(int(r) + r0)
    v_idx = []
    for c in v_cand:
        rs = runs_1d(sub[:, c], gap=2)
        if rs and max(b - a + 1 for a, b in rs) >= min_v:
            v_idx.append(int(c) + c0)

    return dark, gray, clusters(np.array(h_idx, dtype=int)), clusters(np.array(v_idx, dtype=int))


def find_outward_ticks(
    dark: np.ndarray,
    spine: tuple[float, int, int],
    span: tuple[int, int],
    orientation: str,
    side: int = -1,
    band: tuple[int, int] = (4, 16),
    frac: float = 0.55,
    max_width: int = 14,
) -> list[float]:
    """Find outward-pointing ticks along a spine.

    orientation='y': ticks stick horizontally off a vertical spine at `spine`; `span` is (r_min, r_max).
                     `side=-1` looks left of `spine[1]`; `side=+1` looks right of `spine[2]`.
    orientation='x': ticks stick vertically off a horizontal spine at `spine`; `span` is (c_min, c_max).
                     `side=+1` looks below `spine[2]`; `side=-1` looks above `spine[1]`.
    """
    _, s0, s1 = spine
    a, b = span
    b0, b1 = band
    H, W = dark.shape
    if orientation == "y":
        if side < 0:
            c0, c1 = max(0, s0 - b1), max(0, s0 - b0)
        else:
            c0, c1 = min(W, s1 + 1 + b0), min(W, s1 + 1 + b1)
        prof = dark[max(0, a) : min(H, b), c0:c1].mean(axis=1)
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=2)
    elif orientation == "x":
        if side > 0:
            r0, r1 = min(H, s1 + 1 + b0), min(H, s1 + 1 + b1)
        else:
            r0, r1 = max(0, s0 - b1), max(0, s0 - b0)
        prof = dark[r0:r1, max(0, a) : min(W, b)].mean(axis=0)
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=2)
    else:
        raise ValueError(f"Unknown orientation {orientation!r}")
    return [c[0] for c in cs if (c[2] - c[1] + 1) <= max_width]


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 3.5) -> list[float]:
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
            if abs(nm - r_mid) <= max(35.0, 2.5 * thick):
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
    max_search_dx: int = 80,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on an isolated dashed/dotted curve mask via Gaussian-weighted local quadratic fit."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c_lo = max(0, c - max_search_dx)
    c_hi = min(W - 1, c + max_search_dx)
    sub = mask[row_min : row_max + 1, c_lo : c_hi + 1]
    rr, cc = np.nonzero(sub)
    if len(rr) < 12:
        return {"x": float(x_val), "y": None, "found": False, "note": f"too few dash pixels near col {c}"}

    cols = cc.astype(float) + c_lo + 0.5
    rows = rr.astype(float) + row_min + 0.5
    dx = cols - c_f
    sigma = max(20.0, max_search_dx * 0.45)
    w = np.exp(-0.5 * (dx / sigma) ** 2)
    deg = 2 if (np.ptp(dx) > 20 and len(np.unique(cc)) >= 8) else 1
    coef = np.polyfit(dx, rows, deg, w=np.sqrt(w))
    r_est = float(np.polyval(coef, 0.0))
    slope = float(coef[-2]) if len(coef) >= 2 else 0.0

    # Estimate vertical thickness near the fit
    resid = np.abs(rows - np.polyval(coef, dx))
    half_thick = float(np.percentile(resid, 90))
    half_px = max(half_thick, 2.0) + 1.0 + y_axis.residual_px + abs(slope) * (0.5 + x_axis.residual_px)
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
        "note": f"{note_prefix}cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px, dashed_fit_n={len(rr)}",
    }


def make_bar_reading(
    x_val: float,
    col_center: float,
    top_row: float,
    y_axis: Axis,
    upper_err_row: float | None = None,
    lower_err_row: float | None = None,
    mirror_lower: bool = True,
    err_px: float | None = None,
    note: str = "",
) -> dict:
    """Build a calibrated bar-chart reading dict with optional error bar bounds."""
    h_px = err_px if err_px is not None else (1.5 + y_axis.residual_px)
    y_val = float(y_axis.to_value(top_row))
    ya, yb = float(y_axis.to_value(top_row - h_px)), float(y_axis.to_value(top_row + h_px))
    lo, hi = min(ya, yb), max(ya, yb)
    out = {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [float(col_center), float(top_row)],
        "px_lo": float(top_row - h_px),
        "px_hi": float(top_row + h_px),
        "note": note,
    }
    if upper_err_row is not None or lower_err_row is not None:
        u_row = upper_err_row if upper_err_row is not None else top_row
        bar_hi = float(y_axis.to_value(u_row))
        bar_hi_out = float(y_axis.to_value(u_row - h_px))
        if lower_err_row is not None:
            l_row = lower_err_row
            bar_lo = float(y_axis.to_value(l_row))
            bar_lo_out = float(y_axis.to_value(l_row + h_px))
            px_bot = l_row
        elif mirror_lower:
            bar_lo = 2.0 * y_val - bar_hi
            bar_lo_out = 2.0 * y_val - bar_hi_out
            px_bot = top_row + (top_row - u_row)
            out["note"] = (out["note"] + "; " if out["note"] else "") + "lower error bar hidden/unplotted: mirrored from upper"
        else:
            bar_lo = y_val
            bar_lo_out = lo
            px_bot = top_row
        out["bar_lo"] = min(bar_lo, bar_hi)
        out["bar_hi"] = max(bar_lo, bar_hi)
        out["bar_lo_outer"] = min(bar_lo_out, bar_hi_out)
        out["bar_hi_outer"] = max(bar_lo_out, bar_hi_out)
        # Show full error bar span on overlay so visual check verifies error-bar ends too
        out["px_lo"] = float(min(u_row, px_bot))
        out["px_hi"] = float(max(u_row, px_bot))
    return out


def decode_colorbar_rgb(
    img: np.ndarray,
    cb_col_range: tuple[int, int],
    cb_row_top: float,
    cb_row_bot: float,
    val_top: float,
    val_bot: float,
    scale: str = "linear",
    n_samples: int = 1001,
):
    """Sample a vertical colorbar's RGB profile between cb_row_top (val_top) and cb_row_bot (val_bot).

    Returns (lut_vals, lut_rgb) where lut_rgb has shape (n_samples, 3).
    """
    rows = np.linspace(cb_row_top, cb_row_bot, n_samples)
    c0, c1 = cb_col_range
    rgb = np.zeros((n_samples, 3), dtype=float)
    H = img.shape[0]
    for i, r in enumerate(rows):
        ri = int(round(r))
        ri = max(0, min(H - 1, ri))
        rgb[i] = np.median(img[ri, c0:c1, :].astype(float), axis=0)
    if scale == "log10":
        u = np.linspace(math.log10(val_top), math.log10(val_bot), n_samples)
        vals = 10.0 ** u
    else:
        vals = np.linspace(val_top, val_bot, n_samples)
    return vals, rgb


def match_colorbar(rgb_val: np.ndarray, lut_vals: np.ndarray, lut_rgb: np.ndarray) -> tuple[float, float, float]:
    """Match an RGB triplet to the closest colorbar value.

    Returns (best_val, min_rgb_dist, half_width_est) where half_width_est covers all LUT
    entries within min_dist + 4 RGB units (or at least 0.5% of range).
    """
    dists = np.linalg.norm(lut_rgb - rgb_val[None, :], axis=1)
    idx = int(np.argmin(dists))
    dmin = float(dists[idx])
    near = lut_vals[dists <= dmin + 5.0]
    lo = float(near.min())
    hi = float(near.max())
    val = float(lut_vals[idx])
    hw = max(0.5 * (hi - lo), 0.005 * abs(float(lut_vals[0] - lut_vals[-1])))
    return val, dmin, hw
