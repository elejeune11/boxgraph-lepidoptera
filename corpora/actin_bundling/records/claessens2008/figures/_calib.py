"""Shared axis-frame detection, tick calibration, and uncertainty helpers for claessens2008 figures."""
from __future__ import annotations

import math
import numpy as np
from boxgraph.digitize import Axis


def find_frame(img: np.ndarray, gray_thresh: float = 100.0, min_frac: float = 0.45):
    """Detect horizontal and vertical dark frame lines spanning at least min_frac of width/height.
    Returns (h_lines, v_lines), each a list of (start_idx, end_idx, center_coord).
    """
    gray = img.astype(float).mean(axis=2)
    dark = gray < gray_thresh

    def _runs(mask: np.ndarray):
        idx = np.nonzero(mask)[0]
        if len(idx) == 0:
            return []
        splits = np.where(np.diff(idx) > 1)[0] + 1
        return [(int(g[0]), int(g[-1]), float(g.mean() + 0.5)) for g in np.split(idx, splits)]

    return _runs(dark.mean(axis=1) >= min_frac), _runs(dark.mean(axis=0) >= min_frac)


def find_ticks(
    img: np.ndarray,
    border_pos: float,
    span_lo: int,
    span_hi: int,
    axis: str,
    direction: int,
    offset_lo: int = 8,
    offset_hi: int = 24,
    min_depth_frac: float = 0.5,
    min_width: int = 1,
    max_width: int = 18,
    gray_thresh: float = 110.0,
):
    """Detect neutral-dark inward ticks protruding from a frame border.
    axis='x' scans columns along a horizontal border (border_pos is row).
    axis='y' scans rows along a vertical border (border_pos is col).
    direction is +1 (down/right) or -1 (up/left).
    Returns list of (center_px, width_px, depth_px).
    """
    gray = img.astype(float).mean(axis=2)
    spread = img.max(axis=2).astype(float) - img.min(axis=2).astype(float)
    neutral_dark = (gray < gray_thresh) & (spread < 45.0)

    offsets = [int(round(border_pos + direction * d)) for d in range(offset_lo, offset_hi + 1)]
    if axis == "x":
        strip = neutral_dark[offsets, span_lo:span_hi]
        score = strip.mean(axis=0)
        limit = img.shape[0]
    else:
        strip = neutral_dark[span_lo:span_hi, offsets]
        score = strip.mean(axis=1)
        limit = img.shape[1]

    idx = np.nonzero(score >= min_depth_frac)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > 1)[0] + 1
    out = []
    for g in np.split(idx, splits):
        w = len(g)
        if min_width <= w <= max_width:
            c_mid = int(round(g.mean())) + span_lo
            depth = 0
            for d in range(offset_lo, 95):
                p = int(round(border_pos + direction * d))
                if not (0 <= p < limit):
                    break
                val = neutral_dark[p, c_mid] if axis == "x" else neutral_dark[c_mid, p]
                if val:
                    depth = d
                elif d > offset_lo + 4:
                    break
            out.append((float(g.mean() + 0.5 + span_lo), w, depth))
    return out


def value_interval(y_axis: Axis, row_center: float, half_span_px: float) -> tuple[float, float, float]:
    """Return (y, lo, hi) from a continuous row center and pixel half-uncertainty."""
    y = float(y_axis.to_value(row_center))
    a = float(y_axis.to_value(row_center - half_span_px))
    b = float(y_axis.to_value(row_center + half_span_px))
    lo, hi = min(a, b), max(a, b)
    if not (lo <= y <= hi):
        lo, hi = min(lo, y), max(hi, y)
    return y, lo, hi


def sample_curve_at(
    mask: np.ndarray,
    x_axis: Axis,
    y_axis: Axis,
    x_val: float,
    r_lo: int,
    r_hi: int,
    window: int = 2,
    extra_px: float = 1.0,
    note_prefix: str = "",
) -> dict:
    """Read a clean single-curve binary mask at x_val across columns [c-window, c+window]."""
    c_exact = x_axis.to_pixel(x_val)
    c_int = int(round(c_exact))
    H, W = mask.shape
    c0, c1 = max(0, c_int - window), min(W - 1, c_int + window)
    rows_by_col = []
    cols_used = []
    for c in range(c0, c1 + 1):
        rs = np.nonzero(mask[r_lo:r_hi, c])[0]
        if len(rs) > 0:
            rows_by_col.append(rs + r_lo)
            cols_used.append(c)
    if not rows_by_col:
        return {"x": float(x_val), "y": None, "note": f"{note_prefix}no curve pixels at x={x_val:g}"}

    # Take the column closest to c_exact for vertical thickness, and local slope across window
    medians = [float(rs.mean() + 0.5) for rs in rows_by_col]
    if len(cols_used) >= 2:
        slope_rc, intercept_rc = np.polyfit(np.array(cols_used, dtype=float) + 0.5, medians, 1)
        row_c = float(slope_rc * c_exact + intercept_rc)
    else:
        slope_rc = 0.0
        row_c = medians[0]

    # Closest column vertical thickness
    closest_idx = int(np.argmin([abs(c + 0.5 - c_exact) for c in cols_used]))
    rs_mid = rows_by_col[closest_idx]
    thick_px = float(rs_mid.max() - rs_mid.min() + 1)
    # 95% uncertainty half-width in pixels: half stroke thickness + calibration misfit + x-misfit*slope + extra
    cal_y = max(y_axis.residual_px, 0.5)
    cal_x = max(x_axis.residual_px, 0.5)
    half_px = max(0.5 * thick_px, 1.5) + cal_y + abs(slope_rc) * (0.5 + cal_x) + extra_px
    y, lo, hi = value_interval(y_axis, row_c, half_px)
    return {
        "x": float(x_val),
        "y": y,
        "lo": lo,
        "hi": hi,
        "px": [float(c_exact), float(row_c)],
        "px_lo": float(row_c - half_px),
        "px_hi": float(row_c + half_px),
        "note": f"{note_prefix}stroke {thick_px:.0f} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
    }
