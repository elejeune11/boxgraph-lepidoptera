"""Shared fast, image-calibrated figure digitization helpers for winkelman2016."""

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


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 2.5) -> Axis:
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
    gray = image.astype(np.float32).mean(axis=2)
    dark = gray < thresh
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
        rs = runs_1d(sub[r], gap=2)
        if rs and max(b - a + 1 for a, b in rs) >= min_h:
            h_idx.append(int(r) + r0)
    v_idx = []
    for c in v_cand:
        rs = runs_1d(sub[:, c], gap=2)
        if rs and max(b - a + 1 for a, b in rs) >= min_v:
            v_idx.append(int(c) + c0)

    return dark, gray, clusters(np.array(h_idx, dtype=int)), clusters(np.array(v_idx, dtype=int))


def find_ticks(
    dark: np.ndarray,
    spine: tuple[float, int, int],
    span: tuple[int, int],
    orientation: str,
    side: int = -1,
    band: tuple[int, int] = (4, 18),
    frac: float = 0.5,
    max_width: int = 16,
    gap: int = 2,
) -> list[float]:
    """Find tick marks protruding from a spine on `side` (-1 = left/up, +1 = right/down)."""
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
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=gap)
    elif orientation == "x":
        if side > 0:
            r0, r1 = min(H, s1 + 1 + b0), min(H, s1 + 1 + b1)
        else:
            r0, r1 = max(0, s0 - b1), max(0, s0 - b0)
        prof = dark[r0:r1, max(0, a) : min(W, b)].mean(axis=0)
        cs = clusters(np.where(prof >= frac)[0] + max(0, a), gap=gap)
    else:
        raise ValueError(f"Unknown orientation {orientation!r}")
    return [c[0] for c in cs if (c[2] - c[1] + 1) <= max_width]


def match_ticks_regular(
    candidates: list[float], n_expected: int, tol: float = 5.0, min_step: float = 20.0
) -> list[float]:
    """Select the n_expected most evenly spaced tick positions from candidates."""
    c = sorted(candidates)
    if len(c) == n_expected:
        return c
    best = None
    best_err = float("inf")
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = (c[j] - c[i]) / (n_expected - 1)
            if step < min_step:
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
            if abs(nm - r_mid) <= max(40.0, 3.0 * thick):
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


def make_point_reading(
    x_val: float,
    col_center: float,
    row_center: float,
    y_axis: Axis,
    upper_err_row: float | None = None,
    lower_err_row: float | None = None,
    mirror_lower: bool = True,
    err_px: float | None = None,
    note: str = "",
) -> dict:
    """Build a calibrated bar or marker reading dict with optional vertical error bar bounds."""
    h_px = err_px if err_px is not None else (1.5 + y_axis.residual_px)
    y_val = float(y_axis.to_value(row_center))
    ya, yb = float(y_axis.to_value(row_center - h_px)), float(y_axis.to_value(row_center + h_px))
    lo, hi = min(ya, yb), max(ya, yb)
    out = {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [float(col_center), float(row_center)],
        "px_lo": float(row_center - h_px),
        "px_hi": float(row_center + h_px),
        "note": note,
    }
    if upper_err_row is not None or lower_err_row is not None:
        u_row = upper_err_row if upper_err_row is not None else row_center
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
            px_bot = row_center + (row_center - u_row)
            out["note"] = (out["note"] + "; " if out["note"] else "") + "lower error bar hidden/unplotted: mirrored from upper"
        else:
            bar_lo = y_val
            bar_lo_out = lo
            px_bot = row_center
        out["bar_lo"] = min(bar_lo, bar_hi)
        out["bar_hi"] = max(bar_lo, bar_hi)
        out["bar_lo_outer"] = min(bar_lo_out, bar_hi_out)
        out["bar_hi_outer"] = max(bar_lo_out, bar_hi_out)
        out["px_lo"] = float(min(u_row, px_bot))
        out["px_hi"] = float(max(u_row, px_bot))
    return out


def find_vertical_error_bar(
    stem_mask: np.ndarray,
    col_center: float,
    row_center: float,
    row_min: int,
    row_max: int,
    stem_half_w: int = 8,
    max_gap: int = 12,
    cap_half_w: int = 18,
    min_cap_width: int = 10,
    start_up_offset: int = 0,
    start_dn_offset: int = 0,
) -> tuple[float | None, float | None]:
    """Trace a vertical error bar upward and downward from (col_center, row_center).

    Returns (upper_row, lower_row) where each is the center row of the horizontal cap (if found)
    or the extreme end of the contiguous vertical stem (if no cap or cap is small).
    """
    H, W = stem_mask.shape
    c0 = max(0, int(round(col_center)) - stem_half_w)
    c1 = min(W, int(round(col_center)) + stem_half_w + 1)
    strip = stem_mask[:, c0:c1].any(axis=1)

    r_start_up = max(row_min, int(round(row_center)) - start_up_offset)
    r_start_dn = min(row_max, int(round(row_center)) + start_dn_offset)

    # Trace up
    top_r = None
    gap = 0
    for r in range(r_start_up, row_min - 1, -1):
        if strip[r]:
            top_r = r
            gap = 0
        else:
            gap += 1
            if gap > max_gap:
                break

    # Trace down
    bot_r = None
    gap = 0
    for r in range(r_start_dn, row_max + 1):
        if strip[r]:
            bot_r = r
            gap = 0
        else:
            gap += 1
            if gap > max_gap:
                break

    def refine_cap(end_r: int | None, direction: int) -> float | None:
        if end_r is None:
            return None
        # Look within a few pixels of end_r for a wider horizontal cap centered near col_center
        k0 = max(0, int(round(col_center)) - cap_half_w)
        k1 = min(W, int(round(col_center)) + cap_half_w + 1)
        lo_r = max(row_min, end_r - 6)
        hi_r = min(row_max, end_r + 6)
        widths = stem_mask[lo_r : hi_r + 1, k0:k1].sum(axis=1)
        cap_rows = np.where(widths >= min_cap_width)[0] + lo_r
        if len(cap_rows) > 0:
            # Take the cluster closest to the extreme end
            cs = clusters(cap_rows, gap=2)
            best = min(cs, key=lambda t: t[0]) if direction < 0 else max(cs, key=lambda t: t[0])
            return float(best[0])
        return float(end_r) + 0.5

    return refine_cap(top_r, -1), refine_cap(bot_r, +1)
