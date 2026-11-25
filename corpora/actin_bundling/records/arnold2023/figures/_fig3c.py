"""Shared reader for Fig. 3(c) (PDF p. 4 / 128402-4): shear-flow static structure factor <S(k)> vs. k."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from _common import build_axis, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def calibrate_fig3c(image: np.ndarray):
    box = find_box(image)
    c0, r0, c1, r1 = box

    xt = [
        p
        for p, _ in scan_ticks(image, box, "top_in", offset_min=5, offset_max=25, min_len=5, gray_thresh=200)
        if c0 + 50 < p < c1 - 50
    ]
    xa = build_axis("linear", [c0] + xt + [c1], [0.0, 0.5, 1.0, 1.5, 2.0], max_residual_px=2.0)

    yt = [
        p
        for p, _ in match_ticks_both_sides(
            scan_ticks(image, box, "left_in", offset_min=5, offset_max=25, min_len=5, gray_thresh=200),
            scan_ticks(image, box, "right_in", offset_min=5, offset_max=25, min_len=5, gray_thresh=200),
        )
        if r0 + 50 < p < r1 - 50
    ]
    ya = build_axis("linear", [r0] + yt + [r1], [14.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0], max_residual_px=2.0)
    return box, xa, ya


def _clean_curve_mask(image: np.ndarray, box: tuple[float, float, float, float]) -> np.ndarray:
    """Keep the curve+marker network (excluding text labels and frame) and fill only marker-sized holes."""
    c0, r0, c1, r1 = box
    dark = image.astype(np.float64).mean(axis=2) < 140.0
    # Crop strictly inside the box frame (excluding ticks near top/sides, keeping baseline markers just above r1)
    r_lo, r_hi = int(r0) + 28, int(r1) - 4
    c_lo, c_hi = int(c0) + 28, int(c1) - 28
    interior = np.zeros_like(dark)
    interior[r_lo:r_hi, c_lo:c_hi] = dark[r_lo:r_hi, c_lo:c_hi]

    # Erase isolated text glyphs ("t = 50", "t = 30", "t = 20"): keep components wider than 250 px
    lbl, n = ndimage.label(interior)
    slices = ndimage.find_objects(lbl)
    net = np.zeros_like(interior)
    for i, sl in enumerate(slices, start=1):
        if sl is None:
            continue
        if (sl[1].stop - sl[1].start) >= 250:
            net[sl] |= lbl[sl] == i

    # Fill small marker interior white holes (area < 1800 px, w < 70, h < 85) without filling inter-curve loops
    white = ~net
    wlbl, _ = ndimage.label(white)
    wslices = ndimage.find_objects(wlbl)
    solid = net.copy()
    for i, sl in enumerate(wslices, start=1):
        if sl is None:
            continue
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if h <= 85 and w <= 70:
            sub = wlbl[sl] == i
            if int(sub.sum()) <= 1800:
                solid[sl] |= sub
    return solid


def _vertical_clusters(solid: np.ndarray, col: int, r_lo: int, r_hi: int, half_w: int = 2) -> list[tuple[float, float]]:
    """Return sorted (row_center, half_height_px) for distinct vertical runs in solid[:, col-half_w:col+half_w+1]."""
    strip = solid[r_lo:r_hi, col - half_w : col + half_w + 1].any(axis=1)
    rr = np.where(strip)[0] + r_lo
    if len(rr) == 0:
        return []
    gs: list[list[int]] = [[int(rr[0])]]
    for r in rr[1:]:
        if int(r) <= gs[-1][-1] + 6:
            gs[-1].append(int(r))
        else:
            gs.append([int(r)])
    out = []
    for g in gs:
        center = 0.5 * (g[0] + g[-1]) + 0.5
        half_h = max(0.5 * (g[-1] - g[0]), 4.0)
        out.append((center, half_h))
    return out  # sorted top-to-bottom (highest y first)


def _pick_cluster_for_series(series: str, k: float, clusters: list[tuple[float, float]], ya) -> tuple[float, float] | None:
    """Select the vertical cluster belonging to t=50 (diamonds), t=30 (squares), or t=20 (circles) at wavevector k.

    Ordering from top (highest <S(k)>) to bottom at each k (verified directly on the plotted curves):
      k=0.3: top = t=50 (~1.0)
      k=0.4: top = t=50 (~6.4), lower = t=30 (~0.35)
      k=0.5: [t=50 (~12.6), t=30 (~1.3), t=20 (~0.18)]
      k=0.6: [t=50 (~8.9), t=30 (~3.4), t=20 (~0.53)]
      k=0.7: [t=30 (~7.0), t=50 (~5.2), t=20 (~1.2)]
      k=0.8: [t=30 (~5.9), t=50 (~4.1), t=20 (~2.7)]
      k=0.9: [t=20 & t=30 (~3.95, overlapping), t=50 (~2.7)]
      k=1.0: [t=20 (~3.1), t=30 (~2.4), t=50 (~1.5)]
      k=1.1: [t=20 (~2.5), t=30 (~1.45), t=50 (~0.85)]
      k=1.2: [t=20 (~1.7), t=30 (~0.85), t=50 (~0.5)]
      k=1.3: [t=20 (~1.05), t=30 (~0.55), t=50 (~0.25)]
      k=1.4: [t=20 (~0.5), ...]
    """
    if not clusters:
        return None
    # Sort top-to-bottom (descending y)
    cs = sorted(clusters, key=lambda t: t[0])
    n = len(cs)

    if series == "t50":
        if k <= 0.65:
            return cs[0]  # highest
        elif k <= 0.85:
            return cs[1] if n >= 2 else None  # middle (crossed below t=30)
        else:
            return cs[-1]  # lowest (crossed below t=20)
    elif series == "t30":
        if k <= 0.65:
            return cs[1] if n >= 2 else None  # middle
        elif k <= 0.85:
            return cs[0]  # highest
        elif abs(k - 0.9) < 0.03:
            return cs[0]  # overlapping with t=20 at k=0.9
        else:
            return cs[1] if n >= 3 else None  # middle for k >= 1.0
    elif series == "t20":
        if k <= 0.85:
            return cs[-1]  # lowest
        else:
            return cs[0]  # highest for k >= 0.9
    raise ValueError(f"unknown series {series!r}")


def read_series(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    box, xa, ya = calibrate_fig3c(image)
    c0, r0, c1, r1 = box
    solid = _clean_curve_mask(image, box)

    out = []
    for k in at:
        col = xa.to_pixel(k)
        ci = int(round(col))
        clusters = _vertical_clusters(solid, ci, int(r0) + 28, int(r1) - 4, half_w=2)
        picked = _pick_cluster_for_series(series, float(k), clusters, ya)
        if picked is None:
            out.append({"x": float(k), "y": None, "found": False, "note": f"no resolved cluster for {series} at k={k}"})
            continue
        row, half_h = picked
        # Marker center is known to ~quarter of marker height (~10-15 px) unless overlapping at k=0.9
        note = f"Fig. 3(c) {series} marker/curve center ({len(clusters)} curves in col, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)"
        if abs(k - 0.9) < 0.03 and series in ("t20", "t30"):
            note += "; t=20 and t=30 markers overlap at k=0.9"
            eff_half = max(half_h, 25.0)
        else:
            eff_half = min(max(0.45 * half_h, 10.0), 28.0)
        out.append(reading_from_px(xa, ya, col, row, eff_half, x_val=float(k), note=note))
    return out
