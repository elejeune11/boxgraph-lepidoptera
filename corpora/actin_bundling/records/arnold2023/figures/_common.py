"""Shared axis calibration and marker/curve extraction utilities for arnold2023 figures."""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage

from boxgraph import digitize as D


def find_box(image: np.ndarray, gray_thresh: float = 150.0, min_frac: float = 0.45) -> tuple[float, float, float, float]:
    """Detect the main rectangular plot frame (c0, r0, c1, r1) in continuous pixel coordinates."""
    rgb = image.astype(np.float64)
    gray = rgb.mean(axis=2)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    dark = (gray < gray_thresh) & (chroma < 45.0)
    H, W = dark.shape

    def longest_run_lengths(mask_2d: np.ndarray, axis: int, min_len: float) -> np.ndarray:
        n = mask_2d.shape[axis]
        lens = np.zeros(n, dtype=np.int32)
        sums = mask_2d.sum(axis=1 - axis)
        cands = np.where(sums >= min_len)[0]
        for i in cands:
            v = mask_2d[i] if axis == 0 else mask_2d[:, i]
            d = np.diff(np.concatenate(([0], v.view(np.int8), [0])))
            starts = np.where(d == 1)[0]
            ends = np.where(d == -1)[0]
            if len(starts):
                lens[i] = int((ends - starts).max())
        return lens

    r_lens = longest_run_lengths(dark, 0, min_frac * W)
    c_lens = longest_run_lengths(dark, 1, min_frac * H)
    long_r = np.where(r_lens >= min_frac * W)[0]
    long_c = np.where(c_lens >= min_frac * H)[0]

    def group_centers(idxs: np.ndarray, gap: int = 5) -> list[float]:
        if len(idxs) == 0:
            return []
        gs: list[list[int]] = [[int(idxs[0])]]
        for x in idxs[1:]:
            if int(x) <= gs[-1][-1] + gap:
                gs[-1].append(int(x))
            else:
                gs.append([int(x)])
        return [float(np.mean(g)) + 0.5 for g in gs]

    rg = group_centers(long_r)
    cg = group_centers(long_c)
    if len(rg) < 2 or len(cg) < 2:
        raise RuntimeError(f"Could not detect plot box borders: rg={rg}, cg={cg}")
    return cg[0], rg[0], cg[-1], rg[-1]


def scan_ticks(
    image: np.ndarray,
    box: tuple[float, float, float, float],
    side: str,
    offset_min: int = 6,
    offset_max: int = 35,
    min_len: int = 6,
    max_width: int = 20,
    gray_thresh: float = 195.0,
    edge_margin: int = 6,
) -> list[tuple[float, int]]:
    """Find neutral-dark tick marks along one side of the plot box.

    Returns list of (pixel_center, max_length_in_band).
    `side` is one of: bottom_in, bottom_out, top_in, top_out, left_in, left_out, right_in, right_out.
    """
    c0, r0, c1, r1 = box
    rgb = image.astype(np.float64)
    gray = rgb.mean(axis=2)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    dark = (gray < gray_thresh) & (chroma < 45.0)
    H, W = dark.shape

    def clamp_r(a: int, b: int) -> slice:
        return slice(max(0, min(H, a)), max(0, min(H, b)))

    def clamp_c(a: int, b: int) -> slice:
        return slice(max(0, min(W, a)), max(0, min(W, b)))

    if side == "bottom_in":
        rs, cs, axis, base = clamp_r(int(r1) - offset_max, int(r1) - offset_min), clamp_c(int(c0) + edge_margin, int(c1) - edge_margin), 0, max(0, int(c0) + edge_margin)
    elif side == "bottom_out":
        rs, cs, axis, base = clamp_r(int(r1) + offset_min, int(r1) + offset_max), clamp_c(int(c0) - edge_margin, int(c1) + edge_margin), 0, max(0, int(c0) - edge_margin)
    elif side == "top_in":
        rs, cs, axis, base = clamp_r(int(r0) + offset_min, int(r0) + offset_max), clamp_c(int(c0) + edge_margin, int(c1) - edge_margin), 0, max(0, int(c0) + edge_margin)
    elif side == "left_in":
        rs, cs, axis, base = clamp_r(int(r0) + edge_margin, int(r1) - edge_margin), clamp_c(int(c0) + offset_min, int(c0) + offset_max), 1, max(0, int(r0) + edge_margin)
    elif side == "left_out":
        rs, cs, axis, base = clamp_r(int(r0) - edge_margin, int(r1) + edge_margin), clamp_c(int(c0) - offset_max, int(c0) - offset_min), 1, max(0, int(r0) - edge_margin)
    elif side == "right_in":
        rs, cs, axis, base = clamp_r(int(r0) + edge_margin, int(r1) - edge_margin), clamp_c(int(c1) - offset_max, int(c1) - offset_min), 1, max(0, int(r0) + edge_margin)
    else:
        raise ValueError(f"unknown side {side!r}")

    band = dark[rs, cs]
    proj = band.sum(axis=axis)
    hits = np.where(proj >= min_len)[0]
    ticks: list[tuple[float, int]] = []
    if len(hits) == 0:
        return ticks
    g = [int(hits[0])]
    for h in hits[1:]:
        if int(h) <= g[-1] + 4:
            g.append(int(h))
        else:
            if len(g) <= max_width:
                w = proj[g].astype(np.float64)
                pos = float(np.sum(np.array(g, dtype=np.float64) * w) / np.sum(w)) + base + 0.5
                ticks.append((pos, int(w.max())))
            g = [int(h)]
    if len(g) <= max_width:
        w = proj[g].astype(np.float64)
        pos = float(np.sum(np.array(g, dtype=np.float64) * w) / np.sum(w)) + base + 0.5
        ticks.append((pos, int(w.max())))
    return ticks


def match_ticks_both_sides(
    t1: list[tuple[float, int]],
    t2: list[tuple[float, int]],
    tol_px: float = 3.0,
) -> list[tuple[float, int]]:
    """Keep ticks that appear on both opposite axes (averaging their positions)."""
    out: list[tuple[float, int]] = []
    for p1, l1 in t1:
        near = [(p2, l2) for p2, l2 in t2 if abs(p2 - p1) <= tol_px]
        if near:
            p2, l2 = min(near, key=lambda x: abs(x[0] - p1))
            out.append(((p1 + p2) / 2.0, min(l1, l2)))
    return out


def build_axis(scale: str, positions: list[float], values: list[float], max_residual_px: float = 3.5) -> D.Axis:
    """Fit linear or log10 axis calibration and verify low residual."""
    if len(positions) != len(values):
        raise RuntimeError(f"Tick count mismatch: found {len(positions)} positions vs {len(values)} expected values")
    ax = D.Axis(scale, list(zip(positions, values)))
    if ax.residual_px > max_residual_px:
        raise RuntimeError(f"Axis calibration residual too high: {ax.residual_px:.2f} px > {max_residual_px}")
    return ax


def reading_from_px(
    xa: D.Axis,
    ya: D.Axis,
    col: float,
    row: float,
    half_row_px: float,
    x_val: float | None = None,
    note: str = "",
) -> dict:
    """Build a BoxGraph code-mode reading dict from pixel coordinates with honest vertical uncertainty."""
    x = float(x_val) if x_val is not None else float(xa.to_value(col))
    y = float(ya.to_value(row))
    # Combine pixel uncertainty, quantization (0.5 px), and calibration residual in quadrature
    eff_half = math.hypot(max(half_row_px, 1.0), max(ya.residual_px, 0.5)) + 0.5
    v1 = float(ya.to_value(row - eff_half))
    v2 = float(ya.to_value(row + eff_half))
    lo, hi = min(v1, v2), max(v1, v2)
    return {
        "x": x,
        "y": y,
        "lo": lo,
        "hi": hi,
        "px": [float(col), float(row)],
        "px_lo": float(row - eff_half),
        "px_hi": float(row + eff_half),
        "note": note,
    }


def connected_centroids(
    mask: np.ndarray,
    min_area: int = 20,
    max_area: int = 50000,
) -> list[tuple[float, float, int, tuple[int, int, int, int]]]:
    """Return list of (col, row, area, (r_min, c_min, r_max, c_max)) for connected components."""
    lbl, n = ndimage.label(mask)
    if n == 0:
        return []
    slices = ndimage.find_objects(lbl)
    out = []
    for i, sl in enumerate(slices, start=1):
        if sl is None:
            continue
        sub = lbl[sl] == i
        area = int(sub.sum())
        if area < min_area or area > max_area:
            continue
        rr, cc = np.nonzero(sub)
        row = float(rr.mean()) + sl[0].start + 0.5
        col = float(cc.mean()) + sl[1].start + 0.5
        bbox = (sl[0].start, sl[1].start, sl[0].stop, sl[1].stop)
        out.append((col, row, area, bbox))
    return out
