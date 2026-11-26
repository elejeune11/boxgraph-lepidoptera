"""Shared calibration, curve, bar, boxplot, and 3D-projection helpers for durre2018."""

from __future__ import annotations

import numpy as np
from boxgraph.digitize import Axis


def clusters(indices: np.ndarray, gap: int = 3) -> list[tuple[float, int, int]]:
    """Group sorted integer indices into contiguous clusters (center, start, end)."""
    if len(indices) == 0:
        return []
    splits = np.where(np.diff(indices) > gap)[0] + 1
    return [(float(g.mean()) + 0.5, int(g[0]), int(g[-1])) for g in np.split(indices, splits)]


def runs_1d(mask_1d: np.ndarray, gap: int = 1) -> list[tuple[int, int]]:
    """Return (start, end_inclusive) runs of True in a 1D boolean array."""
    idx = np.nonzero(mask_1d)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap + 1)[0] + 1
    return [(int(g[0]), int(g[-1])) for g in np.split(idx, splits)]


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 3.0) -> Axis:
    """Build an Axis and assert its calibration residual is within tolerance."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_residual_px:
        raise RuntimeError(
            f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px: {marks}"
        )
    return ax


def find_box(
    image: np.ndarray,
    thresh: float = 160.0,
    row_frac: float = 0.45,
    col_frac: float = 0.45,
    row_range: tuple[int, int] | None = None,
    col_range: tuple[int, int] | None = None,
):
    """Find closed 4-sided plot box spines: returns (dark, gray, cl, rt, cr, rb).

    Each spine is a tuple (center, start, end) in full-image pixel coordinates.
    """
    gray = image.sum(axis=2, dtype=np.uint16)  # 0..765 (3 * 0..255)
    dark = gray < int(round(thresh * 3.0))
    H, W = dark.shape
    r0, r1 = row_range if row_range is not None else (0, H)
    c0, c1 = col_range if col_range is not None else (0, W)
    sub = dark[r0:r1, c0:c1]

    min_h = row_frac * sub.shape[1]
    min_v = col_frac * sub.shape[0]
    h_idx = []
    for r in np.where(sub.sum(axis=1) >= min_h)[0]:
        rs = runs_1d(sub[r], gap=2)
        if rs and max(b - a + 1 for a, b in rs) >= min_h:
            h_idx.append(int(r) + r0)
    v_idx = []
    for c in np.where(sub.sum(axis=0) >= min_v)[0]:
        rs = runs_1d(sub[:, c], gap=2)
        if rs and max(b - a + 1 for a, b in rs) >= min_v:
            v_idx.append(int(c) + c0)

    hs = clusters(np.array(h_idx, dtype=int), gap=4)
    vs = clusters(np.array(v_idx, dtype=int), gap=4)
    if len(hs) < 2 or len(vs) < 2:
        raise RuntimeError(f"Could not find 4-sided plot box: hs={hs}, vs={vs}")
    return dark, gray, vs[0], hs[0], vs[-1], hs[-1]


def find_inward_ticks(
    dark: np.ndarray,
    cl: tuple[float, int, int],
    rt: tuple[float, int, int],
    cr: tuple[float, int, int],
    rb: tuple[float, int, int],
    side: str,
    band: tuple[int, int] = (4, 18),
    margin: int = 8,
    frac: float = 0.5,
    max_width: int = 24,
) -> list[float]:
    """Find inward-pointing tick centers along 'y_left', 'y_right', 'x_bot', or 'x_top'."""
    b0, b1 = band
    if side == "y_left":
        prof = dark[rt[2] + margin : rb[1] - margin, cl[2] + b0 : cl[2] + b1].mean(axis=1)
        cs = clusters(np.where(prof >= frac)[0] + rt[2] + margin, gap=3)
    elif side == "y_right":
        prof = dark[rt[2] + margin : rb[1] - margin, cr[1] - b1 : cr[1] - b0].mean(axis=1)
        cs = clusters(np.where(prof >= frac)[0] + rt[2] + margin, gap=3)
    elif side == "x_bot":
        prof = dark[rb[1] - b1 : rb[1] - b0, cl[2] + margin : cr[1] - margin].mean(axis=0)
        cs = clusters(np.where(prof >= frac)[0] + cl[2] + margin, gap=3)
    elif side == "x_top":
        prof = dark[rt[2] + b0 : rt[2] + b1, cl[2] + margin : cr[1] - margin].mean(axis=0)
        cs = clusters(np.where(prof >= frac)[0] + cl[2] + margin, gap=3)
    else:
        raise ValueError(f"Unknown side {side!r}")
    return [c[0] for c in cs if (c[2] - c[1] + 1) <= max_width]


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 6.0) -> list[float]:
    """Select the n_expected most evenly spaced positions from candidates."""
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
    search_dx: int = 3,
    extra_px: float = 0.0,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on a boolean curve mask at x_val within [row_min, row_max]."""
    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    c = max(0, min(W - 1, c))
    row_min = max(0, row_min)
    row_max = min(H - 1, row_max)

    def get_runs(col):
        if col < 0 or col >= W:
            return []
        return [(r0 + row_min, r1 + row_min) for r0, r1 in runs_1d(mask[row_min : row_max + 1, col], gap=2)]

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
        return {"x": float(x_val), "y": None, "found": False, "note": f"{note_prefix}no curve pixel at col {c}"}

    if pick_fn is not None:
        seg = pick_fn(rs, used_c)
    else:
        seg = rs[0] if len(rs) == 1 else max(rs, key=lambda t: t[1] - t[0])
    if seg is None:
        return {"x": float(x_val), "y": None, "found": False, "note": f"{note_prefix}no matching run at col {used_c}"}

    r0, r1 = seg
    r_mid = 0.5 * (r0 + r1) + 0.5
    thick = float(r1 - r0 + 1)

    pts = []
    for dc in (-6, -4, -2, 0, 2, 4, 6):
        cc = used_c + dc
        s_rs = get_runs(cc)
        if s_rs:
            near = min(s_rs, key=lambda t: abs(0.5 * (t[0] + t[1]) + 0.5 - r_mid))
            nm = 0.5 * (near[0] + near[1]) + 0.5
            if abs(nm - r_mid) <= max(40.0, 3.0 * thick):
                pts.append((float(cc) + 0.5, nm))
    slope = 0.0
    if len(pts) >= 2:
        xs_p = np.array([p[0] for p in pts])
        ys_p = np.array([p[1] for p in pts])
        if np.ptp(xs_p) > 0:
            slope = float(np.polyfit(xs_p, ys_p, 1)[0])

    r_est = r_mid + slope * (c_f - (used_c + 0.5))
    half_px = max(0.5 * thick, 1.5) + 0.5 + y_axis.residual_px + abs(slope) * (0.5 + x_axis.residual_px) + extra_px
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
    max_gap_px: int = 45,
    pick_fn=None,
    note_prefix: str = "",
) -> dict:
    """Sample y(x) on a dashed curve, interpolating linearly across dash gaps if needed."""
    direct = sample_solid_curve_y(mask, x_axis, y_axis, x_val, row_min, row_max, pick_fn=pick_fn, search_dx=2, note_prefix=note_prefix)
    if direct.get("found", direct.get("y") is not None):
        return direct

    c_f = float(x_axis.to_pixel(x_val))
    c = int(round(c_f))
    H, W = mask.shape
    row_min = max(0, row_min)
    row_max = min(H - 1, row_max)

    def col_mid(col):
        if col < 0 or col >= W:
            return None
        rs = [(r0 + row_min, r1 + row_min) for r0, r1 in runs_1d(mask[row_min : row_max + 1, col], gap=2)]
        if not rs:
            return None
        seg = pick_fn(rs, col) if pick_fn is not None else max(rs, key=lambda t: t[1] - t[0])
        if seg is None:
            return None
        return 0.5 * (seg[0] + seg[1]) + 0.5, float(seg[1] - seg[0] + 1)

    left = None
    for d in range(1, max_gap_px + 1):
        res = col_mid(c - d)
        if res is not None:
            left = (c - d + 0.5, res[0], res[1])
            break
    right = None
    for d in range(1, max_gap_px + 1):
        res = col_mid(c + d)
        if res is not None:
            right = (c + d + 0.5, res[0], res[1])
            break

    if left is not None and right is not None:
        t = (c_f - left[0]) / (right[0] - left[0])
        r_est = left[1] + t * (right[1] - left[1])
        thick = 0.5 * (left[2] + right[2])
        gap_w = right[0] - left[0]
        slope = (right[1] - left[1]) / gap_w
        half_px = max(0.5 * thick, 2.0) + 1.5 + y_axis.residual_px + abs(slope) * (0.15 * gap_w + x_axis.residual_px)
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
            "note": f"{note_prefix}gap_interp({gap_w:.0f}px), cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        }
    return {"x": float(x_val), "y": None, "found": False, "note": f"{note_prefix}no dash within +/-{max_gap_px}px of col {c}"}


def make_point_reading(
    x_val: float,
    col_center: float,
    row_val: float,
    y_axis: Axis,
    half_px: float | None = None,
    bar_row_lo: float | None = None,
    bar_row_hi: float | None = None,
    note: str = "",
) -> dict:
    """Build a calibrated point/bar/boxplot reading with optional whisker/error-bar bounds."""
    h_px = half_px if half_px is not None else (2.0 + y_axis.residual_px)
    y_val = float(y_axis.to_value(row_val))
    ya, yb = float(y_axis.to_value(row_val - h_px)), float(y_axis.to_value(row_val + h_px))
    lo, hi = min(ya, yb), max(ya, yb)
    out = {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [float(col_center), float(row_val)],
        "px_lo": float(row_val - h_px),
        "px_hi": float(row_val + h_px),
        "note": note,
    }
    if bar_row_lo is not None and bar_row_hi is not None:
        v1 = float(y_axis.to_value(bar_row_lo))
        v2 = float(y_axis.to_value(bar_row_hi))
        blo, bhi = min(v1, v2), max(v1, v2)
        b_err = abs(float(y_axis.to_value(row_val - h_px)) - y_val)
        out["bar_lo"] = blo
        out["bar_hi"] = bhi
        out["bar_lo_outer"] = blo - b_err
        out["bar_hi_outer"] = bhi + b_err
    return out
