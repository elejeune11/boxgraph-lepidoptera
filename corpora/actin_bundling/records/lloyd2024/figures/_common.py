"""Shared calibration and pixel-reading helpers for lloyd2024 figure readers."""

from __future__ import annotations

import numpy as np
from boxgraph import digitize as D


def find_runs(indices: np.ndarray, gap: int = 3) -> list[tuple[int, int, float]]:
    """Group sorted 1D integer indices into contiguous runs (start, end, center)."""
    idx = np.asarray(indices, dtype=int)
    if len(idx) == 0:
        return []
    splits = np.split(idx, np.where(np.diff(idx) > gap)[0] + 1)
    return [(int(g[0]), int(g[-1]), float(g.mean() + 0.5)) for g in splits]


def run_centers(indices: np.ndarray, gap: int = 3) -> list[float]:
    return [c for _, _, c in find_runs(indices, gap=gap)]


def find_frame_box(img: np.ndarray, dark_thresh: float = 120.0,
                   min_row_frac: float = 0.55, min_col_frac: float = 0.50,
                   left_margin: int = 0, right_margin: int = 0,
                   top_margin: int = 0, bottom_margin: int = 0,
                   gap: int = 5) -> tuple[float, float, float, float]:
    """Find plot border lines (c0, r0, c1, r1) from dark ink fractions."""
    gray = img.astype(float).mean(axis=2)
    dark = gray < dark_thresh
    H, W = dark.shape
    r_sub = dark[:, left_margin:W - right_margin].mean(axis=1)
    if top_margin:
        r_sub[:top_margin] = 0.0
    if bottom_margin:
        r_sub[H - bottom_margin:] = 0.0
    c_sub = dark[top_margin:H - bottom_margin, :].mean(axis=0)
    if left_margin:
        c_sub[:left_margin] = 0.0
    if right_margin:
        c_sub[W - right_margin:] = 0.0
    rs = run_centers(np.where(r_sub >= min_row_frac)[0], gap=gap)
    cs = run_centers(np.where(c_sub >= min_col_frac)[0], gap=gap)
    if len(rs) < 1 or len(cs) < 1:
        raise ValueError(f"Could not detect axes lines: rows={rs}, cols={cs}")
    return (cs[0], rs[0], cs[-1], rs[-1])


def find_l_axes(img: np.ndarray, dark_thresh: float = 120.0,
                min_row_frac: float = 0.45, min_col_frac: float = 0.45,
                left_margin: int = 40, bottom_margin: int = 0,
                gap: int = 5) -> tuple[float, float]:
    """Find L-shaped left y-axis column c0 and bottom x-axis row r1."""
    gray = img.astype(float).mean(axis=2)
    dark = gray < dark_thresh
    H, W = dark.shape
    c_frac = dark.mean(axis=0)
    c_frac[:left_margin] = 0.0
    r_frac = dark.mean(axis=1)
    if bottom_margin:
        r_frac[H - bottom_margin:] = 0.0
    cs = run_centers(np.where(c_frac >= min_col_frac)[0], gap=gap)
    rs = run_centers(np.where(r_frac >= min_row_frac)[0], gap=gap)
    if not cs or not rs:
        raise ValueError(f"Could not find L-axes: cols={cs}, rows={rs}")
    return (cs[0], rs[-1])


def y_interval(ya: D.Axis, row_center: float, half_px: float,
               extra_val: float = 0.0) -> tuple[float, float, float, float, float]:
    """Convert pixel center and half-thickness to (y, lo, hi, px_lo, px_hi)."""
    h = max(half_px, 1.0) + max(ya.residual_px, 0.5)
    y = ya.to_value(row_center)
    a = ya.to_value(row_center - h)
    b = ya.to_value(row_center + h)
    lo = min(a, b) - extra_val
    hi = max(a, b) + extra_val
    return y, lo, hi, float(row_center - h), float(row_center + h)
