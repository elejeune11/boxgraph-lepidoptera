"""Shared fast axis-calibration and curve/marker digitization helpers for mietke2019.

All calibration comes directly from the rendered panel image `img` passed to `read(img, at)`:
plot-box borders and inward tick marks are detected from neutral dark ink, matched against
the printed tick labels visible on the figure, and least-squares fitted in linear or log10 space.
"""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage


class Axis:
    """Linear or log10 calibration fitted to (pixel, value) tick pairs."""

    def __init__(self, scale: str, marks: list[tuple[float, float]]):
        if scale not in ("linear", "log10"):
            raise ValueError(f"Unsupported axis scale: {scale!r}")
        self.scale = scale
        pts = [(float(p), math.log10(float(v)) if scale == "log10" else float(v)) for p, v in marks]
        if len(pts) < 2:
            raise ValueError(f"Need >= 2 calibration marks, got {len(pts)}")
        ps = np.array([p for p, _ in pts], dtype=float)
        us = np.array([u for _, u in pts], dtype=float)
        self.b, self.a = [float(x) for x in np.polyfit(ps, us, 1)]  # u = a + b * p
        self.residual_px = float(np.max(np.abs(ps - (us - self.a) / self.b)))
        self.marks = list(marks)

    def to_value(self, p: float) -> float:
        u = self.a + self.b * float(p)
        return 10.0 ** u if self.scale == "log10" else u

    def to_pixel(self, v: float) -> float:
        u = math.log10(float(v)) if self.scale == "log10" else float(v)
        return (u - self.a) / self.b


def _clusters(indices: np.ndarray, max_gap: int = 3) -> list[tuple[int, int, float]]:
    if len(indices) == 0:
        return []
    splits = np.where(np.diff(indices) > max_gap)[0] + 1
    out = []
    for grp in np.split(indices, splits):
        out.append((int(grp[0]), int(grp[-1]), float(grp.mean())))
    return out


def neutral_dark_mask(img: np.ndarray, gray_max: float = 115.0, sat_max: float = 35.0) -> np.ndarray:
    arr = img.astype(float)
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    return (gray < gray_max) & (sat < sat_max)


def find_plot_box(img: np.ndarray, min_frac: float = 0.45) -> tuple[float, float, float, float]:
    """Detect rectangular plot box (c0, r0, c1, r1) from neutral dark axis borders."""
    black = neutral_dark_mask(img)
    hc = _clusters(np.where(black.mean(axis=1) > min_frac)[0])
    vc = _clusters(np.where(black.mean(axis=0) > min_frac)[0])
    if len(hc) < 2 or len(vc) < 2:
        raise RuntimeError(f"Could not detect rectangular plot box (H={len(hc)}, V={len(vc)})")
    r0, r1 = hc[0][2], hc[-1][2]
    c0, c1 = vc[0][2], vc[-1][2]
    return (c0, r0, c1, r1)


def detect_interior_ticks(
    img: np.ndarray,
    box: tuple[float, float, float, float],
    side: str,
    end_margin: int = 22,
    band_lo: int = 11,
    band_hi: int = 23,
    min_fill: float = 0.45,
) -> list[float]:
    """Detect centers of inward neutral-dark tick marks along one side of the plot box."""
    c0, r0, c1, r1 = box
    black = neutral_dark_mask(img)
    H, W = black.shape
    if side == "bottom":
        rs = slice(max(0, int(round(r1)) - band_hi), max(0, int(round(r1)) - band_lo))
        cs = slice(int(round(c0)) + end_margin, int(round(c1)) - end_margin)
        prof = black[rs, cs].mean(axis=0)
        offset = int(round(c0)) + end_margin
    elif side == "top":
        rs = slice(min(H, int(round(r0)) + band_lo), min(H, int(round(r0)) + band_hi))
        cs = slice(int(round(c0)) + end_margin, int(round(c1)) - end_margin)
        prof = black[rs, cs].mean(axis=0)
        offset = int(round(c0)) + end_margin
    elif side == "left":
        rs = slice(int(round(r0)) + end_margin, int(round(r1)) - end_margin)
        cs = slice(min(W, int(round(c0)) + band_lo), min(W, int(round(c0)) + band_hi))
        prof = black[rs, cs].mean(axis=1)
        offset = int(round(r0)) + end_margin
    elif side == "right":
        rs = slice(int(round(r0)) + end_margin, int(round(r1)) - end_margin)
        cs = slice(max(0, int(round(c1)) - band_hi), max(0, int(round(c1)) - band_lo))
        prof = black[rs, cs].mean(axis=1)
        offset = int(round(r0)) + end_margin
    else:
        raise ValueError(f"Unknown side: {side!r}")
    return [m + offset for _, _, m in _clusters(np.where(prof > min_fill)[0])]


def match_ticks_linear(
    detected: list[float],
    expected_values: list[float],
    vmin: float,
    vmax: float,
    p_vmin: float,
    p_vmax: float,
    include_endpoints: bool = True,
    tol_frac: float = 0.03,
) -> Axis:
    """Match detected tick pixels to expected printed values along a linear axis.

    `p_vmin` and `p_vmax` give the box border pixels corresponding to `vmin` and `vmax`
    (for a vertical y-axis, `p_vmin = r1` at bottom and `p_vmax = r0` at top).
    """
    span_p = p_vmax - p_vmin
    span_v = vmax - vmin
    tol_px = abs(span_p) * tol_frac
    marks: list[tuple[float, float]] = []
    for val in expected_values:
        frac = (val - vmin) / span_v
        pred = p_vmin + frac * span_p
        if detected:
            best = min(detected, key=lambda p: abs(p - pred))
            if abs(best - pred) <= tol_px:
                marks.append((best, val))
                continue
        if include_endpoints and (abs(val - vmin) < 1e-12 or abs(val - vmax) < 1e-12):
            marks.append((p_vmin if abs(val - vmin) < 1e-12 else p_vmax, val))
    ax = Axis("linear", marks)
    if ax.residual_px > 1.5:
        raise RuntimeError(f"Linear axis calibration residual too large: {ax.residual_px:.2f} px (marks={marks})")
    return ax


def calibrate_shifted_linear_ticks(
    detected: list[float],
    expected_interior_values: list[float],
    ascending_pixel_means_ascending_value: bool = False,
) -> Axis:
    """Calibrate an axis whose border endpoints may be offset from the outer box (e.g. dual y-axis).

    Matches evenly-spaced detected interior ticks directly to `expected_interior_values`.
    """
    if len(detected) < len(expected_interior_values):
        raise RuntimeError(f"Expected >= {len(expected_interior_values)} ticks, detected {detected}")
    det = sorted(detected)
    vals = sorted(expected_interior_values, reverse=not ascending_pixel_means_ascending_value)
    # Find the contiguous or best-fitting evenly spaced subset of length len(vals)
    k = len(vals)
    if len(det) == k:
        best_marks = list(zip(det, vals))
    else:
        from itertools import combinations
        best_ax = None
        best_marks = None
        for comb in combinations(det, k):
            ax = Axis("linear", list(zip(comb, vals)))
            if best_ax is None or ax.residual_px < best_ax.residual_px:
                best_ax = ax
                best_marks = list(zip(comb, vals))
    ax = Axis("linear", best_marks)
    if ax.residual_px > 1.5:
        raise RuntimeError(f"Shifted linear axis calibration residual too large: {ax.residual_px:.2f} px")
    return ax


def match_ticks_log10(
    detected: list[float],
    expected_values: list[float],
    vmin: float,
    vmax: float,
    p_vmin: float,
    p_vmax: float,
    include_endpoints: bool = True,
    tol_frac: float = 0.03,
) -> Axis:
    """Match detected tick pixels to expected printed values along a log10 axis."""
    lmin, lmax = math.log10(vmin), math.log10(vmax)
    span_p = p_vmax - p_vmin
    span_l = lmax - lmin
    tol_px = abs(span_p) * tol_frac
    marks: list[tuple[float, float]] = []
    for val in expected_values:
        frac = (math.log10(val) - lmin) / span_l
        pred = p_vmin + frac * span_p
        if detected:
            best = min(detected, key=lambda p: abs(p - pred))
            if abs(best - pred) <= tol_px:
                marks.append((best, val))
                continue
        if include_endpoints and (abs(math.log10(val) - lmin) < 1e-12 or abs(math.log10(val) - lmax) < 1e-12):
            marks.append((p_vmin if abs(math.log10(val) - lmin) < 1e-12 else p_vmax, val))
    ax = Axis("log10", marks)
    if ax.residual_px > 1.5:
        raise RuntimeError(f"Log10 axis calibration residual too large: {ax.residual_px:.2f} px (marks={marks})")
    return ax


def rgb_mask(img: np.ndarray, rgb: tuple[float, float, float], tol: float) -> np.ndarray:
    target = np.array(rgb, dtype=float)
    d = np.sqrt(((img.astype(float) - target) ** 2).sum(axis=2))
    return d <= tol


def box_interior_mask(shape: tuple[int, int], box: tuple[float, float, float, float], inset: int = 8) -> np.ndarray:
    c0, r0, c1, r1 = box
    H, W = shape
    m = np.zeros((H, W), dtype=bool)
    m[max(0, int(round(r0)) + inset):min(H, int(round(r1)) - inset),
      max(0, int(round(c0)) + inset):min(W, int(round(c1)) - inset)] = True
    return m


def _runs_1d(arr: np.ndarray, gap: int = 2) -> list[tuple[float, float]]:
    """Return (center, half_width) for each True run in a 1D boolean array."""
    idx = np.where(arr)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > gap + 1)[0] + 1
    out = []
    for grp in np.split(idx, splits):
        lo, hi = float(grp[0]), float(grp[-1])
        # Pixel center coordinates: pixel i spans [i, i+1), center i + 0.5
        out.append(((lo + hi + 1.0) / 2.0, (hi - lo + 1.0) / 2.0))
    return out


def read_vertical_curve(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    at: list[float],
    box: tuple[float, float, float, float],
    search: int = 3,
    dashed: bool = False,
    dash_window: int = 40,
    min_run_px: int = 3,
    branch_selector=None,
) -> list[dict]:
    """Read y(x) at each x in `at` from a boolean curve mask."""
    c0, r0, c1, r1 = box
    H, W = mask.shape
    x_mis = max(0.5, x_axis.residual_px)
    y_mis = max(0.5, y_axis.residual_px)
    out = []

    for xv in at:
        pc = x_axis.to_pixel(xv)
        col = int(round(pc))
        if dashed:
            # Collect all valid runs in [col - dash_window, col + dash_window] and fit local line
            cs_list, rs_list, hs_list = [], [], []
            for c in range(max(int(c0) + 6, col - dash_window), min(int(c1) - 6, col + dash_window) + 1):
                runs = [r for r in _runs_1d(mask[:, c]) if 2.0 * r[1] >= min_run_px]
                if not runs:
                    continue
                if branch_selector is not None and len(runs) > 1:
                    r_best = branch_selector(xv, runs, y_axis)
                else:
                    r_best = max(runs, key=lambda r: r[1])
                cs_list.append(c + 0.5)
                rs_list.append(r_best[0])
                hs_list.append(r_best[1])
            if len(cs_list) < 4:
                out.append({"x": float(xv), "y": None, "note": f"no dashed curve pixels near x={xv:g}"})
                continue
            cs_arr = np.array(cs_list)
            rs_arr = np.array(rs_list)
            # Linear local regression around pc
            slope, intercept = np.polyfit(cs_arr - pc, rs_arr, 1)
            center = float(intercept)
            fit_res = float(np.percentile(np.abs(rs_arr - (intercept + slope * (cs_arr - pc))), 90))
            half = float(np.median(hs_list))
            h_px = half + 0.5 + fit_res + abs(float(slope)) * x_mis + y_mis
            note = f"dashed local fit over {len(cs_list)} cols (res={fit_res:.1f} px, cal={y_axis.residual_px:.2f} px)"
        else:
            runs = []
            used = col
            for dc in sorted(range(-search, search + 1), key=abs):
                c = col + dc
                if 0 <= c < W:
                    cand = [r for r in _runs_1d(mask[:, c]) if 2.0 * r[1] >= min_run_px]
                    if cand:
                        runs = cand
                        used = c
                        break
            if not runs:
                out.append({"x": float(xv), "y": None, "note": f"no curve pixels near x={xv:g}"})
                continue
            if branch_selector is not None and len(runs) > 1:
                center, half = branch_selector(xv, runs, y_axis)
            elif len(runs) > 1:
                # Pick widest (thickest) solid run
                center, half = max(runs, key=lambda r: r[1])
            else:
                center, half = runs[0]
            # Estimate local slope over +/- 4 columns
            slope = 0.0
            left_pts, right_pts = [], []
            for d in (2, 3, 4):
                if used - d >= 0:
                    lr = _runs_1d(mask[:, used - d])
                    if lr:
                        left_pts.append((used - d + 0.5, min(lr, key=lambda r: abs(r[0] - center))[0]))
                if used + d < W:
                    rr = _runs_1d(mask[:, used + d])
                    if rr:
                        right_pts.append((used + d + 0.5, min(rr, key=lambda r: abs(r[0] - center))[0]))
            if left_pts and right_pts:
                lc_x, lc_y = np.mean([p[0] for p in left_pts]), np.mean([p[1] for p in left_pts])
                rc_x, rc_y = np.mean([p[0] for p in right_pts]), np.mean([p[1] for p in right_pts])
                slope = float((rc_y - lc_y) / (rc_x - lc_x))
            offset = abs(pc - (used + 0.5))
            h_px = half + 0.5 + abs(slope) * (offset + x_mis) + y_mis
            note = f"half={half:.1f} px, cal={y_axis.residual_px:.2f} px"

        yv = y_axis.to_value(center)
        a = y_axis.to_value(center - h_px)
        b = y_axis.to_value(center + h_px)
        lo, hi = min(a, b), max(a, b)
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(pc), float(center)],
            "px_lo": float(center - h_px),
            "px_hi": float(center + h_px),
            "note": note,
        })
    return out


def read_horizontal_curve(
    mask: np.ndarray,
    alpha_axis: Axis,
    L_axis: Axis,
    at_alpha: list[float],
    box: tuple[float, float, float, float],
    search: int = 3,
    dashed: bool = False,
    dash_window: int = 45,
    min_run_px: int = 3,
    branch_selector=None,
) -> list[dict]:
    """Read horizontal coordinate L(alpha) at each vertical-axis value alpha in `at_alpha`.

    Returns x = alpha, y = L, px = [col_center, row_alpha] (no px_lo/px_hi because overlay_px
    draws px_span vertically rather than horizontally).
    """
    c0, r0, c1, r1 = box
    H, W = mask.shape
    a_mis = max(0.5, alpha_axis.residual_px)
    L_mis = max(0.5, L_axis.residual_px)
    out = []

    for av in at_alpha:
        pr = alpha_axis.to_pixel(av)
        row = int(round(pr))
        if dashed:
            rs_list, cs_list, hs_list = [], [], []
            for r in range(max(int(r0) + 6, row - dash_window), min(int(r1) - 6, row + dash_window) + 1):
                runs = [x for x in _runs_1d(mask[r, :]) if 2.0 * x[1] >= min_run_px]
                if not runs:
                    continue
                if branch_selector is not None and len(runs) > 1:
                    c_best = branch_selector(av, runs, L_axis)
                else:
                    c_best = max(runs, key=lambda x: x[1])
                rs_list.append(r + 0.5)
                cs_list.append(c_best[0])
                hs_list.append(c_best[1])
            if len(rs_list) < 4:
                out.append({"x": float(av), "y": None, "note": f"no dashed curve pixels near alpha={av:g}"})
                continue
            rs_arr = np.array(rs_list)
            cs_arr = np.array(cs_list)
            # Robust median filter against outlier dash fragments / markers
            med_c = float(np.median(cs_arr))
            keep = np.abs(cs_arr - med_c) < 80.0
            rs_arr, cs_arr = rs_arr[keep], cs_arr[keep]
            if len(rs_arr) < 4:
                out.append({"x": float(av), "y": None, "note": f"insufficient dashed curve pixels near alpha={av:g}"})
                continue
            slope, intercept = np.polyfit(rs_arr - pr, cs_arr, 1)
            center = float(intercept)
            fit_res = float(np.percentile(np.abs(cs_arr - (intercept + slope * (rs_arr - pr))), 90))
            half = float(np.median(hs_list))
            h_px = half + 0.5 + fit_res + abs(float(slope)) * a_mis + L_mis
            note = f"dashed horizontal fit over {len(rs_arr)} rows (res={fit_res:.1f} px)"
        else:
            runs = []
            used = row
            for dr in sorted(range(-search, search + 1), key=abs):
                r = row + dr
                if 0 <= r < H:
                    cand = [x for x in _runs_1d(mask[r, :]) if 2.0 * x[1] >= min_run_px]
                    if cand:
                        runs = cand
                        used = r
                        break
            if not runs:
                out.append({"x": float(av), "y": None, "note": f"no horizontal curve pixels near alpha={av:g}"})
                continue
            if branch_selector is not None and len(runs) > 1:
                center, half = branch_selector(av, runs, L_axis)
            elif len(runs) > 1:
                center, half = max(runs, key=lambda x: x[1])
            else:
                center, half = runs[0]
            # Local slope dc/dr
            slope = 0.0
            up_pts, dn_pts = [], []
            for d in (2, 3, 4):
                if used - d >= 0:
                    ur = _runs_1d(mask[used - d, :])
                    if ur:
                        up_pts.append((used - d + 0.5, min(ur, key=lambda x: abs(x[0] - center))[0]))
                if used + d < H:
                    dr_ = _runs_1d(mask[used + d, :])
                    if dr_:
                        dn_pts.append((used + d + 0.5, min(dr_, key=lambda x: abs(x[0] - center))[0]))
            if up_pts and dn_pts:
                u_r, u_c = np.mean([p[0] for p in up_pts]), np.mean([p[1] for p in up_pts])
                d_r, d_c = np.mean([p[0] for p in dn_pts]), np.mean([p[1] for p in dn_pts])
                slope = float((d_c - u_c) / (d_r - u_r))
            offset = abs(pr - (used + 0.5))
            h_px = half + 0.5 + abs(slope) * (offset + a_mis) + L_mis
            note = f"half={half:.1f} px, cal={L_axis.residual_px:.2f} px"

        Lv = L_axis.to_value(center)
        a = L_axis.to_value(center - h_px)
        b = L_axis.to_value(center + h_px)
        lo, hi = min(a, b), max(a, b)
        out.append({
            "x": float(av),
            "y": float(Lv),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(center), float(pr)],
            "note": note,
        })
    return out


def read_horizontal_region_right_edge(
    region_mask: np.ndarray,
    alpha_axis: Axis,
    L_axis: Axis,
    at_alpha: list[float],
    box: tuple[float, float, float, float],
    min_width: int = 8,
) -> list[dict]:
    """Read the rightmost edge L(alpha) of a shaded region at each vertical-axis alpha."""
    c0, r0, c1, r1 = box
    H, W = region_mask.shape
    a_mis = max(0.5, alpha_axis.residual_px)
    L_mis = max(0.5, L_axis.residual_px)
    out = []

    def right_edge_at_row(r: int) -> float | None:
        if r < int(r0) + 5 or r > int(r1) - 5:
            return None
        runs = _runs_1d(region_mask[r, :], gap=5)
        wide = [x for x in runs if 2.0 * x[1] >= min_width]
        if not wide:
            return None
        # Rightmost wide run's right edge
        c_ctr, hw = max(wide, key=lambda x: x[0] + x[1])
        return c_ctr + hw

    for av in at_alpha:
        pr = alpha_axis.to_pixel(av)
        row = int(round(pr))
        edges = []
        for dr in range(-4, 5):
            e = right_edge_at_row(row + dr)
            if e is not None:
                edges.append((row + dr + 0.5, e))
        if len(edges) < 3:
            out.append({"x": float(av), "y": None, "note": f"no shaded region edge near alpha={av:g}"})
            continue
        rs_arr = np.array([p[0] for p in edges])
        es_arr = np.array([p[1] for p in edges])
        slope, intercept = np.polyfit(rs_arr - pr, es_arr, 1)
        center = float(intercept)
        fit_res = float(np.max(np.abs(es_arr - (intercept + slope * (rs_arr - pr)))))
        # Edge anti-aliasing + slope + calibration uncertainty
        h_px = 3.0 + fit_res + abs(float(slope)) * a_mis + L_mis
        Lv = L_axis.to_value(center)
        a = L_axis.to_value(center - h_px)
        b = L_axis.to_value(center + h_px)
        out.append({
            "x": float(av),
            "y": float(Lv),
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [float(center), float(pr)],
            "note": f"shaded right edge (fit_res={fit_res:.1f} px)",
        })
    return out


def read_markers_xo(
    img: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    box: tuple[float, float, float, float],
    kind: str,
    exclude_rects: list[tuple[float, float, float, float]],
    area_range: tuple[int, int] = (500, 1600),
    size_range: tuple[int, int] = (32, 62),
) -> list[dict]:
    """Detect black 'x' (filled center) or 'o' (hollow center) plot markers."""
    c0, r0, c1, r1 = box
    arr = img.astype(float)
    gray = arr.mean(axis=2)
    sat = arr.max(axis=2) - arr.min(axis=2)
    dark = (gray < 130.0) & (sat < 40.0)
    mask = dark & box_interior_mask(dark.shape, box, inset=25)
    for ec0, er0, ec1, er1 in exclude_rects:
        mask[int(er0):int(er1), int(ec0):int(ec1)] = False

    lbl, n = ndimage.label(mask)
    y_mis = max(0.5, y_axis.residual_px)
    found = []
    for k in range(1, n + 1):
        rr, cc = np.where(lbl == k)
        area = len(rr)
        h = int(rr.max() - rr.min() + 1)
        w = int(cc.max() - cc.min() + 1)
        if not (area_range[0] <= area <= area_range[1] and size_range[0] <= h <= size_range[1] and size_range[0] <= w <= size_range[1]):
            continue
        cy, cx = float(rr.mean()) + 0.5, float(cc.mean()) + 0.5
        ir, ic = int(round(cy)), int(round(cx))
        cdark = float((gray[ir - 3:ir + 4, ic - 3:ic + 4] < 130.0).mean())
        mk = "x" if cdark > 0.4 else ("o" if cdark < 0.2 else "other")
        if mk != kind:
            continue
        xv = x_axis.to_value(cx)
        yv = y_axis.to_value(cy)
        # Honest 95% interval for marker centroid: 1/4 of marker height + 1 px + axis misfit
        h_px = 0.25 * h + 1.0 + y_mis
        a = y_axis.to_value(cy - h_px)
        b = y_axis.to_value(cy + h_px)
        found.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [float(cx), float(cy)],
            "px_lo": float(cy - h_px),
            "px_hi": float(cy + h_px),
            "note": f"marker '{kind}' {w}x{h} px (area={area}, cal_x={x_axis.residual_px:.2f} px, cal_y={y_axis.residual_px:.2f} px)",
        })
    found.sort(key=lambda r: r["x"])
    return found
