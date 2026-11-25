"""Shared calibration, curve-sampling, and marker-detection helpers for You et al. (PNAS 2020) figures."""

from __future__ import annotations

import math
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis


def clusters(indices: np.ndarray, gap: int = 3) -> list[tuple[float, int, int]]:
    """Group sorted 1D integer indices into (center, start, end_inclusive)."""
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


def find_box(
    image: np.ndarray,
    thresh: float = 140.0,
    chroma_max: float = 45.0,
    row_frac: float = 0.45,
    col_frac: float = 0.45,
):
    """Detect the rectangular plot frame (h_top, h_bot, v_left, v_right)."""
    gray = image.astype(float).mean(axis=2)
    chroma = image.max(axis=2).astype(float) - image.min(axis=2).astype(float)
    dark = (gray < thresh) & (chroma < chroma_max)
    h_lines = clusters(np.where(dark.mean(axis=1) > row_frac)[0])
    v_lines = clusters(np.where(dark.mean(axis=0) > col_frac)[0])
    if len(h_lines) < 2 or len(v_lines) < 2:
        raise RuntimeError(
            f"Could not detect rectangular plot box: h_lines={h_lines}, v_lines={v_lines}"
        )
    return dark, gray, h_lines[0], h_lines[-1], v_lines[0], v_lines[-1]


def find_inward_ticks(
    dark: np.ndarray,
    h_top: tuple[float, int, int],
    h_bot: tuple[float, int, int],
    v_left: tuple[float, int, int],
    v_right: tuple[float, int, int],
    band: tuple[int, int] = (4, 20),
    max_width: int = 16,
    frac: float = 0.50,
    margin: int = 3,
):
    """Detect inward-pointing ticks along the bottom and left frame borders."""
    b0, b1 = band
    c0, c1 = v_left[2] + margin, v_right[1] - margin
    r0, r1 = h_top[2] + margin, h_bot[1] - margin

    xt_raw = clusters(
        np.where(dark[h_bot[1] - b1 : h_bot[1] - b0, c0:c1].mean(axis=0) > frac)[0] + c0
    )
    yt_raw = clusters(
        np.where(dark[r0:r1, v_left[2] + b0 : v_left[2] + b1].mean(axis=1) > frac)[0] + r0
    )
    xt = [c[0] for c in xt_raw if (c[2] - c[1] + 1) <= max_width]
    yt = [r[0] for r in yt_raw if (r[2] - r[1] + 1) <= max_width]
    return xt, yt


def match_ticks_regular(candidates: list[float], n_expected: int, tol: float = 3.0) -> list[float]:
    """Select the n_expected most evenly spaced tick positions from candidates."""
    c = sorted(candidates)
    if len(c) == n_expected:
        return c
    if len(c) < n_expected:
        raise RuntimeError(f"Expected {n_expected} regular ticks, only found {len(c)}: {c}")
    best = None
    best_err = float("inf")
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            step = (c[j] - c[i]) / (n_expected - 1)
            if step < 20.0:
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


def verify_axis(scale: str, marks: list[tuple[float, float]], max_residual_px: float = 1.5) -> Axis:
    """Construct an Axis and assert its linear calibration residual is small."""
    ax = Axis(scale, marks)
    if ax.residual_px > max_residual_px:
        raise RuntimeError(
            f"Axis calibration residual {ax.residual_px:.2f} px exceeds {max_residual_px} px: {marks}"
        )
    return ax


def rgb_mask(
    image: np.ndarray,
    target_rgb: tuple[float, float, float],
    tol: float = 75.0,
    min_chroma: float = 30.0,
) -> np.ndarray:
    """Boolean mask of pixels within Euclidean RGB distance tol of target_rgb."""
    arr = image.astype(float)
    t = np.array(target_rgb, dtype=float)
    d = np.sqrt(np.sum((arr - t) ** 2, axis=2))
    mask = d <= tol
    if min_chroma > 0.0:
        chroma = arr.max(axis=2) - arr.min(axis=2)
        mask &= chroma >= min_chroma
    return mask


def mask_interior(
    mask: np.ndarray,
    h_top: tuple[float, int, int],
    h_bot: tuple[float, int, int],
    v_left: tuple[float, int, int],
    v_right: tuple[float, int, int],
    inset: int = 5,
) -> np.ndarray:
    """Restrict a boolean mask strictly to the interior of the plot box."""
    out = np.zeros_like(mask, dtype=bool)
    r0 = h_top[2] + inset
    r1 = h_bot[1] - inset
    c0 = v_left[2] + inset
    c1 = v_right[1] - inset
    out[r0 : r1 + 1, c0 : c1 + 1] = mask[r0 : r1 + 1, c0 : c1 + 1]
    return out


def sample_curve_y(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    row_min: int,
    row_max: int,
    pick_fn=None,
    search_dc: int = 4,
    extra_px: float = 0.0,
    note_prefix: str = "",
) -> dict:
    """Sample y(x_val) on a curve mask within rows [row_min, row_max]."""
    c_f = float(x_axis.to_pixel(x_val))
    c0 = int(round(c_f))
    H, W = mask.shape
    c0 = max(0, min(W - 1, c0))

    rs = []
    c_used = c0
    for dc in sorted(range(-search_dc, search_dc + 1), key=abs):
        cc = c0 + dc
        if 0 <= cc < W:
            sub = mask[row_min : row_max + 1, cc]
            cand = [(a + row_min, b + row_min) for a, b in runs_1d(sub, gap=2)]
            if cand:
                rs = cand
                c_used = cc
                break

    if not rs:
        return {"x": float(x_val), "y": None, "note": f"{note_prefix}no curve pixels near x={x_val:g}"}

    if pick_fn is not None:
        chosen = pick_fn(rs)
        if chosen is None:
            return {"x": float(x_val), "y": None, "note": f"{note_prefix}selector rejected runs at x={x_val:g}"}
        r_start, r_end = chosen
    else:
        # default: widest contiguous run (or single run)
        r_start, r_end = max(rs, key=lambda ab: ab[1] - ab[0])

    r_ctr = 0.5 * (r_start + r_end) + 0.5
    half_stroke = 0.5 * (r_end - r_start + 1)

    # Estimate local slope d(row)/d(col) across nearby columns
    slope = 0.0
    nb = []
    for k in (-3, -2, 2, 3):
        cc = c_used + k
        if 0 <= cc < W:
            sub = mask[row_min : row_max + 1, cc]
            cand = [(a + row_min, b + row_min) for a, b in runs_1d(sub, gap=2)]
            if cand:
                near = min(cand, key=lambda ab: abs(0.5 * (ab[0] + ab[1]) + 0.5 - r_ctr))
                nb.append((cc, 0.5 * (near[0] + near[1]) + 0.5))
    if len(nb) >= 2:
        xs_nb = np.array([p[0] for p in nb], dtype=float)
        ys_nb = np.array([p[1] for p in nb], dtype=float)
        if xs_nb.max() > xs_nb.min():
            slope = float(np.polyfit(xs_nb, ys_nb, 1)[0])

    # Interpolate center slightly if c_used != c_f
    r_eval = r_ctr + slope * (c_f - (c_used + 0.5))
    x_misfit = max(x_axis.residual_px, 0.5)
    y_misfit = max(y_axis.residual_px, 0.5)
    h_px = half_stroke + 0.5 + abs(slope) * (abs(c_f - (c_used + 0.5)) + x_misfit) + y_misfit + extra_px

    y_val = float(y_axis.to_value(r_eval))
    v1 = float(y_axis.to_value(r_eval - h_px))
    v2 = float(y_axis.to_value(r_eval + h_px))
    lo, hi = min(v1, v2), max(v1, v2)
    return {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [float(c_f), float(r_eval)],
        "px_lo": float(r_eval - h_px),
        "px_hi": float(r_eval + h_px),
        "note": f"{note_prefix}stroke_half={half_stroke:.1f}px, slope={slope:.2f}, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
    }


def make_marker_reading(
    x_val: float,
    r_ctr: float,
    c_ctr: float,
    x_axis: Axis,
    y_axis: Axis,
    half_px: float = 2.0,
    note: str = "",
) -> dict:
    """Build a reading dict from a detected marker centroid (c_ctr, r_ctr)."""
    y_misfit = max(y_axis.residual_px, 0.5)
    h_px = half_px + 0.5 + y_misfit
    y_val = float(y_axis.to_value(r_ctr))
    v1 = float(y_axis.to_value(r_ctr - h_px))
    v2 = float(y_axis.to_value(r_ctr + h_px))
    lo, hi = min(v1, v2), max(v1, v2)
    return {
        "x": float(x_val),
        "y": y_val,
        "lo": lo,
        "hi": hi,
        "px": [float(c_ctr), float(r_ctr)],
        "px_lo": float(r_ctr - h_px),
        "px_hi": float(r_ctr + h_px),
        "note": note,
    }


def connected_blobs(mask: np.ndarray, min_area: int = 40, max_area: int = 20000):
    """Label connected components in mask, returning list of dicts sorted by cx."""
    lab, n = ndimage.label(mask)
    blobs = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        area = int(len(rr))
        if area < min_area or area > max_area:
            continue
        blobs.append({
            "cx": float(cc.mean()) + 0.5,
            "cy": float(rr.mean()) + 0.5,
            "r0": int(rr.min()),
            "r1": int(rr.max()),
            "c0": int(cc.min()),
            "c1": int(cc.max()),
            "area": area,
        })
    return sorted(blobs, key=lambda b: b["cx"])
