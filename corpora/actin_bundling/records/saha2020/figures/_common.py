"""Shared fast calibration and curve-reading helpers for saha2020 figures.

All axis calibrations and curve extractions operate directly on the RGB `image`
array passed by BoxGraph, with no PDF re-rendering and no hardcoded pixel
coordinates.
"""

from __future__ import annotations

import math
import numpy as np


def _group_indices(indices: np.ndarray, max_gap: int = 3) -> list[tuple[int, int, float]]:
    """Group sorted 1D integer indices into (start, end, center) runs."""
    if len(indices) == 0:
        return []
    out: list[tuple[int, int, float]] = []
    s = prev = int(indices[0])
    for v in indices[1:]:
        iv = int(v)
        if iv <= prev + max_gap:
            prev = iv
        else:
            out.append((s, prev, 0.5 * (s + prev) + 0.5))
            s = prev = iv
    out.append((s, prev, 0.5 * (s + prev) + 0.5))
    return out


def find_frame(image: np.ndarray, thr: float = 190.0, min_frac: float = 0.45,
               neutral_only: bool = False, max_sat: float = 45.0):
    """Locate the main rectangular plot frame (left, top, right, bottom spines).

    Returns ((c0a, c0b, c0), (r0a, r0b, r0), (c1a, c1b, c1), (r1a, r1b, r1)).
    """
    arr = image.astype(np.float32)
    gray = arr.mean(axis=2)
    mask = gray < thr
    if neutral_only:
        sat = arr.max(axis=2) - arr.min(axis=2)
        mask &= sat < max_sat
    r_idx = np.where(mask.mean(axis=1) > min_frac)[0]
    c_idx = np.where(mask.mean(axis=0) > min_frac)[0]
    rgs = _group_indices(r_idx, max_gap=4)
    cgs = _group_indices(c_idx, max_gap=4)
    if len(rgs) < 2 or len(cgs) < 2:
        raise RuntimeError(f"Could not locate plot frame: rows={rgs}, cols={cgs}")
    # Choose the widest horizontal pair and tallest vertical pair (excluding tiny edge artifacts)
    r_top, r_bot = rgs[0], rgs[-1]
    c_left, c_right = cgs[0], cgs[-1]
    return c_left, r_top, c_right, r_bot


def slice_ticks(mask_2d: np.ndarray, fixed_pos: int, lo: int, hi: int,
                is_x: bool, max_width: int = 25) -> list[float]:
    """Find tick center coordinates along a 1D slice parallel to a spine."""
    line = mask_2d[fixed_pos, lo:hi] if is_x else mask_2d[lo:hi, fixed_pos]
    idx = np.where(line)[0]
    gs = _group_indices(idx, max_gap=3)
    return [g[2] + lo for g in gs if (g[1] - g[0] + 1) <= max_width]


def inward_ticks(mask_2d: np.ndarray, spine_inner: int, direction: int,
                 lo: int, hi: int, is_x: bool, min_len: int = 8,
                 max_len: int = 75, max_width: int = 22) -> list[tuple[float, int]]:
    """Measure inward tick run lengths from a spine; returns [(center, length)]."""
    H, W = mask_2d.shape
    span = max_len + 6
    if is_x:
        rs = [spine_inner + direction * d for d in range(1, span + 1)]
        rs = [r for r in rs if 0 <= r < H]
        strip = mask_2d[rs, lo:hi]  # (D, N)
    else:
        cs = [spine_inner + direction * d for d in range(1, span + 1)]
        cs = [c for c in cs if 0 <= c < W]
        strip = mask_2d[lo:hi, cs].T  # (D, N)
    # Count consecutive True starting at depth 0 or 1
    D, N = strip.shape
    lengths = np.zeros(N, dtype=int)
    for j in range(N):
        col = strip[:, j]
        k = 0
        while k < min(3, D) and not col[k]:
            k += 1
        n = 0
        while k + n < D and col[k + n]:
            n += 1
        lengths[j] = n
    valid = np.where((lengths >= min_len) & (lengths <= max_len))[0]
    gs = _group_indices(valid, max_gap=3)
    out = []
    for s, e, ctr in gs:
        if (e - s + 1) <= max_width:
            out.append((ctr + lo, int(lengths[s:e + 1].max())))
    return out


class LinearAxis:
    """Linear or log10 axis calibration from (pixel, value) pairs."""

    def __init__(self, marks: list[tuple[float, float]], scale: str = "linear"):
        self.scale = scale
        ps = np.array([float(p) for p, _ in marks], dtype=float)
        vs = np.array([math.log10(float(v)) if scale == "log10" else float(v) for _, v in marks], dtype=float)
        self.b, self.a = np.polyfit(ps, vs, 1)  # u = a + b * p
        preds = (vs - self.a) / self.b
        self.residual_px = float(np.max(np.abs(ps - preds))) if len(ps) > 2 else 0.5
        self.misfit_px = max(self.residual_px, 0.75)

    def to_pixel(self, val: float) -> float:
        u = math.log10(val) if self.scale == "log10" else val
        return float((u - self.a) / self.b)

    def to_value(self, px: float) -> float:
        u = self.a + self.b * px
        return float(10.0 ** u) if self.scale == "log10" else float(u)


def column_runs(mask: np.ndarray, col: int, r_lo: int, r_hi: int,
                max_gap: int = 2) -> list[tuple[int, int, float, float]]:
    """Return [(start, end, center, half_height)] of True runs in column `col`."""
    if not (0 <= col < mask.shape[1]):
        return []
    r0 = max(0, r_lo)
    r1 = min(mask.shape[0], r_hi)
    idx = np.where(mask[r0:r1, col])[0] + r0
    gs = _group_indices(idx, max_gap=max_gap)
    return [(s, e, ctr, 0.5 * (e - s + 1)) for s, e, ctr in gs]


def estimate_slope(mask: np.ndarray, col: int, center: float, r_lo: int, r_hi: int,
                   delta: int = 3) -> float:
    """Estimate local pixel slope dr/dc near (col, center)."""
    left = column_runs(mask, col - delta, r_lo, r_hi)
    right = column_runs(mask, col + delta, r_lo, r_hi)
    if not left or not right:
        return 0.0
    lc = min(left, key=lambda t: abs(t[2] - center))[2]
    rc = min(right, key=lambda t: abs(t[2] - center))[2]
    return float((rc - lc) / (2.0 * delta))


def make_reading(xv: float, pc: float, center: float, half_px: float, slope: float,
                 xa: LinearAxis, ya: LinearAxis, extra_px: float = 0.5, note: str = "") -> dict:
    """Build a BoxGraph code-mode reading dict with honest 95% pixel/value interval."""
    col_offset = abs(pc - (math.floor(pc) + 0.5))
    h_px = half_px + extra_px + abs(slope) * (col_offset + xa.misfit_px) + ya.misfit_px
    y = ya.to_value(center)
    v1 = ya.to_value(center - h_px)
    v2 = ya.to_value(center + h_px)
    lo, hi = (v1, v2) if v1 <= v2 else (v2, v1)
    return {
        "x": float(xv),
        "y": float(y),
        "lo": float(lo),
        "hi": float(hi),
        "px": [float(pc), float(center)],
        "px_lo": float(center - h_px),
        "px_hi": float(center + h_px),
        "note": note,
    }


def read_branch_from_mask(mask: np.ndarray, xa: LinearAxis, ya: LinearAxis,
                          at: list[float], r_lo: int, r_hi: int,
                          pick: str = "single", search: int = 3,
                          min_run_px: int = 1, max_gap: int = 2,
                          extra_px: float = 0.5, note_prefix: str = "") -> list[dict]:
    """Sample a curve mask at x positions `at`.

    `pick` selects which vertical run to take when multiple branches cross a column:
      - 'single' / 'thickest': widest/thickest run
      - 'top': smallest row (highest y value)
      - 'bottom': largest row (lowest y value)
      - 'second_from_top': second smallest row (e.g. right-upper hyperbola below top branch)
      - 'second_from_bottom': second largest row
      - 'nearest:<row>': run closest to given row
    """
    out = []
    for xv in at:
        pc = xa.to_pixel(xv)
        col0 = int(round(pc))
        runs = []
        used = col0
        for dc in sorted(range(-search, search + 1), key=abs):
            cand = [t for t in column_runs(mask, col0 + dc, r_lo, r_hi, max_gap=max_gap)
                    if (t[1] - t[0] + 1) >= min_run_px]
            if cand:
                runs = cand
                used = col0 + dc
                break
        if not runs:
            out.append({"x": float(xv), "y": None, "note": f"no mask pixels near x={xv:g}"})
            continue
        if pick == "top":
            chosen = runs[0]
        elif pick == "bottom":
            chosen = runs[-1]
        elif pick == "second_from_top":
            chosen = runs[1] if len(runs) >= 2 else runs[-1]
        elif pick == "second_from_bottom":
            chosen = runs[-2] if len(runs) >= 2 else runs[0]
        elif pick == "thickest" or pick == "single":
            chosen = max(runs, key=lambda t: t[3])
        elif pick.startswith("nearest:"):
            target = float(pick.split(":", 1)[1])
            chosen = min(runs, key=lambda t: abs(t[2] - target))
        else:
            chosen = runs[0]
        _, _, ctr, half = chosen
        slope = estimate_slope(mask, used, ctr, r_lo, r_hi)
        note = f"{note_prefix}; cal_misfit=({xa.residual_px:.2f},{ya.residual_px:.2f})px"
        out.append(make_reading(xv, pc, ctr, half, slope, xa, ya, extra_px=extra_px, note=note.strip("; ")))
    return out
