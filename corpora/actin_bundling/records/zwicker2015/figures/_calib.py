"""Shared axis calibration and curve/marker extraction utilities for zwicker2015 figures.

All panel crops are rendered at 1200 dpi from corpus/actin_bundling/2015_Zwicker.pdf.
Every axis is calibrated directly from the black plot-box spines and inward tick marks
detected in the rendered crop image.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy import ndimage

from boxgraph import digitize as D


def find_plot_box(img: np.ndarray, black_thresh: int = 70) -> dict:
    """Detect the main rectangular plot frame (4 black spines) in a 1200 dpi crop.

    Returns a dict with:
      c_left, r_top, c_right, r_bot: continuous spine center coordinates
      left_in, top_in, right_in, bot_in: integer interior boundary pixels
      h_interior_lines: list of continuous row centers of full-width black interior lines (e.g. y=0)
    """
    black = img.max(axis=2) < black_thresh
    H, W = black.shape

    def long_runs(mask2d: np.ndarray, min_len: int):
        out = []
        for i in range(mask2d.shape[0]):
            v = mask2d[i]
            if int(v.sum()) < min_len:
                continue
            idx = np.nonzero(v)[0]
            splits = np.where(np.diff(idx) > 4)[0] + 1
            best = max(np.split(idx, splits), key=lambda s: s[-1] - s[0] + 1)
            L = int(best[-1] - best[0] + 1)
            if L >= min_len:
                out.append((i, L, int(best[0]), int(best[-1])))
        return out

    def cluster(lines):
        gs = []
        for it in lines:
            if not gs or it[0] > gs[-1][-1][0] + 3:
                gs.append([it])
            else:
                gs[-1].append(it)
        res = []
        for g in gs:
            lo, hi = g[0][0], g[-1][0]
            s = int(np.median([x[2] for x in g]))
            e = int(np.median([x[3] for x in g]))
            res.append((lo, hi, (lo + hi) / 2.0 + 0.5, e - s + 1, s, e))
        return res

    hg = cluster(long_runs(black, int(0.30 * W)))
    vg = cluster(long_runs(black.T, int(0.30 * H)))
    if len(hg) < 2 or len(vg) < 2:
        raise RuntimeError("Could not detect rectangular plot spines")

    max_hl = max(x[3] for x in hg)
    hg_main = [x for x in hg if x[3] >= 0.85 * max_hl]
    top_sp, bot_sp = hg_main[0], hg_main[-1]
    left_sp = min(vg, key=lambda x: abs(x[2] - top_sp[4]))
    right_sp = min(vg, key=lambda x: abs(x[2] - top_sp[5]))

    return {
        "c_left": left_sp[2],
        "r_top": top_sp[2],
        "c_right": right_sp[2],
        "r_bot": bot_sp[2],
        "left_lo": left_sp[0],
        "left_hi": left_sp[1],
        "right_lo": right_sp[0],
        "right_hi": right_sp[1],
        "top_lo": top_sp[0],
        "top_hi": top_sp[1],
        "bot_lo": bot_sp[0],
        "bot_hi": bot_sp[1],
        "left_in": left_sp[1] + 2,
        "right_in": right_sp[0] - 2,
        "top_in": top_sp[1] + 2,
        "bot_in": bot_sp[0] - 2,
        "h_interior_lines": [x[2] for x in hg_main[1:-1]],
    }


def find_zero_line(img: np.ndarray, box: dict, black_thresh: int = 70) -> float:
    """Find a full-width black horizontal line inside the plot box (e.g. y = 0)."""
    black = img.max(axis=2) < black_thresh
    c0, c1 = box["left_in"] + 10, box["right_in"] - 10
    r0, r1 = box["top_in"] + 20, box["bot_in"] - 20
    sub = black[r0:r1, c0:c1]
    frac = sub.mean(axis=1)
    rows = np.nonzero(frac > 0.70)[0]
    if len(rows) == 0:
        raise RuntimeError("Could not find horizontal zero line inside plot box")
    splits = np.where(np.diff(rows) > 3)[0] + 1
    best = max(np.split(rows, splits), key=lambda s: float(frac[s].sum()))
    return float(r0 + (best[0] + best[-1]) / 2.0 + 0.5)


def find_inward_ticks(
    img: np.ndarray,
    box: dict,
    side: str,
    min_len: int = 8,
    max_len: int = 55,
    max_width: int = 14,
    black_thresh: int = 70,
) -> list[tuple[float, int]]:
    """Detect inward black tick marks along 'bottom', 'top', 'left', or 'right'.

    Returns sorted list of (continuous_pixel_position, tick_length_px).
    """
    black = img.max(axis=2) < black_thresh
    H, W = black.shape
    depth = max_len + 10
    c0, c1 = box["left_in"], box["right_in"]
    r0, r1 = box["top_in"], box["bot_in"]

    if side == "bottom":
        strip = black[box["bot_lo"] - 1 : max(0, box["bot_lo"] - 1 - depth) : -1, c0:c1]
        offset = c0
    elif side == "top":
        strip = black[box["top_hi"] + 1 : min(H, box["top_hi"] + 1 + depth), c0:c1]
        offset = c0
    elif side == "left":
        strip = black[r0:r1, box["left_hi"] + 1 : min(W, box["left_hi"] + 1 + depth)].T
        offset = r0
    elif side == "right":
        strip = black[r0:r1, box["right_lo"] - 1 : max(0, box["right_lo"] - 1 - depth) : -1].T
        offset = r0
    else:
        raise ValueError(f"Unknown side {side!r}")

    L = np.zeros(strip.shape[1], dtype=int)
    for j in range(strip.shape[1]):
        col = strip[:, j]
        k = 0
        while k < 2 and k < len(col) and not col[k]:
            k += 1
        n = 0
        while k + n < len(col) and col[k + n]:
            n += 1
        L[j] = n

    pos = np.where((L >= min_len) & (L <= max_len))[0]
    ticks: list[tuple[float, int]] = []
    if len(pos) > 0:
        for s in np.split(pos, np.where(np.diff(pos) > 2)[0] + 1):
            if len(s) <= max_width:
                ticks.append((float(offset + (s[0] + s[-1]) / 2.0 + 0.5), int(L[s].max())))
    return ticks


def filter_ticks_by_len(
    ticks: list[tuple[float, int]], lo: int = 30, hi: int = 48
) -> list[float]:
    return [p for p, L in ticks if lo <= L <= hi]


def match_ticks_regular(
    candidates: list[float], n_expected: int, tol_px: float = 3.0
) -> list[float]:
    """Select the n_expected evenly spaced tick positions from candidates."""
    c = sorted(candidates)
    if len(c) == n_expected:
        return c
    if len(c) < n_expected:
        raise RuntimeError(f"Expected {n_expected} ticks, found only {len(c)}: {c}")
    best = None
    best_err = 1e9
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = c[j] - c[i]
            if step < 30.0:
                continue
            pts = [c[i], c[j]]
            ok = True
            for k in range(2, n_expected):
                target = c[i] + k * step
                near = min(c, key=lambda x: abs(x - target))
                if abs(near - target) <= tol_px:
                    pts.append(near)
                else:
                    ok = False
                    break
            if ok:
                ks = np.arange(n_expected, dtype=float)
                slope, intercept = np.polyfit(ks, pts, 1)
                err = float(np.max(np.abs(np.array(pts) - (intercept + slope * ks))))
                if err < best_err:
                    best_err = err
                    best = pts
    if best is None:
        raise RuntimeError(f"Could not find {n_expected} regular ticks among {c}")
    return best


def fit_axis(scale: str, marks: list[tuple[float, float]], max_res: float = 1.5) -> D.Axis:
    ax = D.Axis(scale, marks)
    if ax.residual_px > max_res:
        raise RuntimeError(
            f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_res} px; marks={marks}"
        )
    return ax


def verify_log_minor_ticks(
    ax: D.Axis,
    all_ticks: list[tuple[float, int]],
    decade_values: list[float],
    tol_px: float = 2.0,
) -> D.Axis:
    """Refine a log10 axis calibration by matching all detected major and minor ticks."""
    marks: list[tuple[float, float]] = []
    pos_list = [p for p, _ in all_ticks]
    for dec in decade_values:
        for m in range(1, 10):
            val = m * dec
            pred = ax.to_pixel(val)
            near = [p for p in pos_list if abs(p - pred) <= tol_px]
            if near:
                best = min(near, key=lambda p: abs(p - pred))
                marks.append((best, val))
    if len(marks) >= 3:
        refined = D.Axis("log10", marks)
        if refined.residual_px <= 1.5:
            return refined
    return ax


# --------------------------------------------------------------------------
# Color masks tuned for APS / Matplotlib palette in Zwicker et al. (2015)
# --------------------------------------------------------------------------


def mask_blue(img: np.ndarray) -> np.ndarray:
    """Dark/medium blue curve ink (matplotlib C0-like: R<80, G~90-170, B>140)."""
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    return (b > 130) & (b - r > 55) & (r < 120)


def mask_orange(img: np.ndarray) -> np.ndarray:
    """Orange / light-orange curve or marker fill (R>200, G~80-210, B<145, R-B>75)."""
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    return (r > 190) & (r - b > 70) & (g > 60) & (g < 220)


def mask_green(img: np.ndarray) -> np.ndarray:
    """Green curve or marker fill (G>110, G-R>35, G-B>35)."""
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    return (g > 100) & (g - r > 35) & (g - b > 35)


def mask_red(img: np.ndarray) -> np.ndarray:
    """Red dashed line (Fig. 9b: R>180, G<70, B<70)."""
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    b = img[:, :, 2].astype(float)
    return (r > 170) & (g < 90) & (b < 90) & (r - g > 100)


def mask_black(img: np.ndarray, thresh: int = 75) -> np.ndarray:
    """Pure black ink (R,G,B all below thresh)."""
    return img.max(axis=2) < thresh


# --------------------------------------------------------------------------
# Curve and marker reading helpers
# --------------------------------------------------------------------------


def _runs_1d(indices: np.ndarray, max_gap: int = 1) -> list[tuple[int, int]]:
    if len(indices) == 0:
        return []
    splits = np.where(np.diff(indices) > max_gap + 1)[0] + 1
    return [(int(s[0]), int(s[-1])) for s in np.split(indices, splits)]


def _col_runs(mask: np.ndarray, col: int, min_run: int = 2, max_run: int = 60) -> list[tuple[float, float]]:
    if not (0 <= col < mask.shape[1]):
        return []
    rows = np.nonzero(mask[:, col])[0]
    out = []
    for s, e in _runs_1d(rows):
        h = e - s + 1
        if min_run <= h <= max_run:
            out.append(((s + e) / 2.0 + 0.5, h / 2.0))
    return out


def interior_mask(mask: np.ndarray, box: dict, inset: int = 8) -> np.ndarray:
    out = np.zeros_like(mask, dtype=bool)
    r0 = max(0, box["top_in"] + inset)
    r1 = min(mask.shape[0], box["bot_in"] - inset)
    c0 = max(0, box["left_in"] + inset)
    c1 = min(mask.shape[1], box["right_in"] - inset)
    out[r0:r1, c0:c1] = mask[r0:r1, c0:c1]
    return out


def make_reading(
    xv: float,
    pc: float,
    row_center: float,
    half_px: float,
    slope_px: float,
    x_axis: D.Axis,
    y_axis: D.Axis,
    col_used: float | None = None,
    extra_px: float = 0.0,
    note: str = "",
) -> dict:
    offset = 0.0 if col_used is None else abs(pc - (col_used + 0.5))
    mx = max(x_axis.residual_px, 0.5)
    my = max(y_axis.residual_px, 0.5)
    h = half_px + 0.5 + abs(slope_px) * (offset + mx) + my + extra_px
    y = y_axis.to_value(row_center)
    a = y_axis.to_value(row_center - h)
    b = y_axis.to_value(row_center + h)
    return {
        "x": float(xv),
        "y": float(y),
        "lo": float(min(a, b)),
        "hi": float(max(a, b)),
        "px": [float(pc), float(row_center)],
        "px_lo": float(row_center - h),
        "px_hi": float(row_center + h),
        "note": note,
    }


def read_curve_at(
    mask: np.ndarray,
    x_axis: D.Axis,
    y_axis: D.Axis,
    at: list[float],
    search: int = 3,
    min_run: int = 2,
    max_run: int = 80,
    select: str | Callable = "single",
    extra_px: float = 0.0,
    slope_step: int = 3,
) -> list[dict]:
    """Read a single-valued curve from a binary mask at requested x values.

    `select`:
      'single': pick the run closest to median of nearby columns (or top if tied)
      'top': smallest row coordinate (highest y value)
      'bottom': largest row coordinate (lowest y value)
      callable(xv, runs, y_axis) -> (center, half)
    """
    out = []
    for xv in at:
        pc = x_axis.to_pixel(xv)
        col = int(math.floor(pc))
        runs = []
        used = col
        for dc in sorted(range(-search, search + 1), key=abs):
            r = _col_runs(mask, col + dc, min_run=min_run, max_run=max_run)
            if r:
                runs = r
                used = col + dc
                break
        if not runs:
            out.append({"x": float(xv), "y": None, "note": f"no curve pixels near x={xv:g}"})
            continue
        if callable(select):
            center, half = select(xv, runs, y_axis)
        elif select == "top":
            center, half = min(runs, key=lambda t: t[0])
        elif select == "bottom":
            center, half = max(runs, key=lambda t: t[0])
        else:
            if len(runs) == 1:
                center, half = runs[0]
            else:
                neigh = [
                    c
                    for k in (-4, -3, -2, 2, 3, 4)
                    for c, _ in _col_runs(mask, used + k, min_run=min_run, max_run=max_run)
                ]
                ref = float(np.median(neigh)) if neigh else runs[0][0]
                center, half = min(runs, key=lambda t: abs(t[0] - ref))
        # Estimate local slope d(row)/d(col)
        left = _col_runs(mask, used - slope_step, min_run=min_run, max_run=max_run)
        right = _col_runs(mask, used + slope_step, min_run=min_run, max_run=max_run)
        slope = 0.0
        if left and right:
            lc = min(left, key=lambda t: abs(t[0] - center))[0]
            rc = min(right, key=lambda t: abs(t[0] - center))[0]
            slope = (rc - lc) / (2.0 * slope_step)
        out.append(
            make_reading(
                xv,
                pc,
                center,
                half,
                slope,
                x_axis,
                y_axis,
                col_used=used,
                extra_px=extra_px,
            )
        )
    return out


def read_dotted_curve_at(
    mask: np.ndarray,
    x_axis: D.Axis,
    y_axis: D.Axis,
    at: list[float],
    min_area: int = 12,
    max_area: int = 450,
    max_aspect: float = 3.0,
    extra_px: float = 1.5,
) -> list[dict]:
    """Read a dotted curve by detecting individual dots and interpolating cleanly at `at`."""
    lab, n = ndimage.label(mask)
    dots = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        area = len(rr)
        if not (min_area <= area <= max_area):
            continue
        h = int(rr.max() - rr.min() + 1)
        w = int(cc.max() - cc.min() + 1)
        if max(h, w) / max(1, min(h, w)) > max_aspect:
            continue
        dots.append((float(cc.mean() + 0.5), float(rr.mean() + 0.5), h / 2.0))
    dots.sort(key=lambda d: d[0])
    if len(dots) < 2:
        raise RuntimeError(f"Found only {len(dots)} dots on dotted curve")
    xs = np.array([d[0] for d in dots])
    ys = np.array([d[1] for d in dots])
    hs = np.array([d[2] for d in dots])

    out = []
    for xv in at:
        pc = x_axis.to_pixel(xv)
        # Check if a dot sits close to pc
        j = int(np.argmin(np.abs(xs - pc)))
        if abs(xs[j] - pc) <= 8.0:
            # Local slope from neighbor dots
            j0 = max(0, j - 1)
            j1 = min(len(dots) - 1, j + 1)
            slope = (ys[j1] - ys[j0]) / max(1.0, xs[j1] - xs[j0])
            row_est = ys[j] + slope * (pc - xs[j])
            out.append(
                make_reading(
                    xv,
                    pc,
                    row_est,
                    float(hs[j]),
                    slope,
                    x_axis,
                    y_axis,
                    col_used=None,
                    extra_px=extra_px,
                    note="from dot centroid",
                )
            )
            continue
        # Otherwise bracket between consecutive dots
        idx = int(np.searchsorted(xs, pc))
        if idx == 0 or idx >= len(dots):
            out.append({"x": float(xv), "y": None, "note": f"x={xv:g} outside dotted curve range"})
            continue
        x0, y0, h0 = xs[idx - 1], ys[idx - 1], hs[idx - 1]
        x1, y1, h1 = xs[idx], ys[idx], hs[idx]
        t = (pc - x0) / (x1 - x0)
        slope = (y1 - y0) / (x1 - x0)
        row_est = y0 + t * (y1 - y0)
        # Curvature safety margin if 3 dots available
        curv_extra = extra_px
        if idx >= 2 and idx < len(dots):
            # Compare quadratic vs linear interpolation
            k0 = max(0, idx - 2)
            k1 = min(len(dots), k0 + 3)
            if k1 - k0 == 3:
                poly = np.polyfit(xs[k0:k1], ys[k0:k1], 2)
                row_quad = float(np.polyval(poly, pc))
                curv_extra += abs(row_quad - row_est)
                row_est = row_quad
        out.append(
            make_reading(
                xv,
                pc,
                row_est,
                float(max(h0, h1)),
                slope,
                x_axis,
                y_axis,
                col_used=None,
                extra_px=curv_extra,
                note="interpolated between adjacent dots",
            )
        )
    return out


def read_filled_markers(
    mask: np.ndarray,
    x_axis: D.Axis,
    y_axis: D.Axis,
    min_area: int = 120,
    max_area: int = 3500,
    shape: str | None = None,
    split_touching: bool = False,
    extra_px: float = 1.0,
) -> list[dict]:
    """Detect filled markers (circles or squares) and return sorted readings.

    `shape`: None (all), 'circle' (fill_ratio < 0.88), or 'square' (fill_ratio >= 0.88).
    """
    work = mask.copy()
    if split_touching:
        work = ndimage.binary_erosion(work, structure=np.ones((3, 3), dtype=bool), iterations=1)
    lab, n = ndimage.label(work)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        area = len(rr)
        if not (min_area <= area <= max_area):
            continue
        h = int(rr.max() - rr.min() + 1)
        w = int(cc.max() - cc.min() + 1)
        aspect = max(h, w) / max(1, min(h, w))
        if aspect > 1.45:
            continue
        fill = area / float(h * w)
        if shape == "circle" and not (0.60 <= fill < 0.87):
            continue
        if shape == "square" and fill < 0.87:
            continue
        cx = float(cc.mean() + 0.5)
        cy = float(rr.mean() + 0.5)
        xv = x_axis.to_value(cx)
        my = max(y_axis.residual_px, 0.5)
        err_px = 1.5 + my + extra_px
        yv = y_axis.to_value(cy)
        a = y_axis.to_value(cy - err_px)
        b = y_axis.to_value(cy + err_px)
        pts.append(
            {
                "x": float(xv),
                "y": float(yv),
                "lo": float(min(a, b)),
                "hi": float(max(a, b)),
                "px": [cx, cy],
                "px_lo": float(cy - err_px),
                "px_hi": float(cy + err_px),
                "note": f"marker area={area} px, fill={fill:.2f}",
            }
        )
    pts.sort(key=lambda d: d["x"])
    return pts
