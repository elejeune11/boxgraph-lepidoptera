"""Shared calibration, histogram marker reader, and log-normal fit curve reader for Fig. 2D (PDF p. 4, 1200 dpi)."""

from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (
    LinearCal,
    detect_ticks_along_spine,
    find_black_spines,
    read_colored_curve_at,
    select_major_ticks,
    split_runs,
)


def _calibrate_fig2d(image: np.ndarray):
    im = image.astype(float)
    H, W, _ = im.shape
    gray = im.mean(axis=2)
    dark = gray < 105.0
    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: 0, 1, 2, 3, 4, 5 um
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=60, min_len=15, coord_min=vl[1], coord_max=W - 10)
    xt = select_major_ticks(xt_raw, expected=6, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], max_res_px=2.5)

    # Y axis: 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=60, min_len=15, coord_min=20, coord_max=hr[0] - 5)
    yt = select_major_ticks(yt_raw, expected=7, frac=0.72)
    y_cal = LinearCal(yt, [0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0], max_res_px=2.5)
    return hr, vl, x_cal, y_cal


def _series_mask(im: np.ndarray, series: str) -> np.ndarray:
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    if series == "005":
        return (gray > 50) & (gray < 125) & (sat < 22)
    elif series == "0125":
        return (R > 185) & (R - G > 110) & (R - B > 110)
    elif series == "025":
        return (B > 155) & (B - R > 110) & (B - G > 40)
    elif series == "050":
        return (G > 125) & (G - R > 60) & (G - B > 25)
    else:
        raise ValueError(series)


def read_fig2d_hist(image: np.ndarray, series: str, bin_centers: list[float]) -> list[dict]:
    """Read discrete histogram markers at exact 0.2 um bin centers in Fig. 2D."""
    im = image.astype(float)
    H, W, _ = im.shape
    hr, vl, x_cal, y_cal = _calibrate_fig2d(image)

    mask = _series_mask(im, series)
    # Mask out legend box (x > 2.8, y > 0.38)
    mask[: int(y_cal.to_px(0.37)), int(x_cal.to_px(2.75)) :] = False
    mask[hr[0] - 2 :, :] = False
    mask[:, : vl[1] + 4] = False

    # Separate markers from thin fit curves (~10 px wide) via 12x12 opening (preserves partially occluded markers!)
    k = 12
    markers = ndimage.binary_opening(mask, structure=np.ones((k, k), dtype=bool))

    # Also compute marker centroid offset correction for triangles:
    # Up-triangle centroid is at 1/3 height from base; plotted data point in OriginLab is the bounding-box center (or centroid)!
    # Wait: OriginLab plots up/down triangle symbols centered at the midpoint of their bounding box (center of circumscribed square/circle),
    # whereas centroid of an equilateral triangle is shifted by ~1/6 height toward the flat base. Let's use bounding-box vertical midpoint
    # for clean unoccluded triangles/squares/circles, or verify on overlay!
    out = []
    half_bin_px = int(0.08 * abs(x_cal.to_px(1.0) - x_cal.to_px(0.0)))  # ~55 px around bin center
    for xv in bin_centers:
        pc = x_cal.to_px(xv)
        c0 = max(vl[1] + 5, int(round(pc)) - half_bin_px)
        c1 = min(W - 5, int(round(pc)) + half_bin_px)
        sub = markers[20 : hr[0] - 2, c0:c1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 40:
            # Fallback for heavily occluded gray square behind red circle at x=1.2: check raw mask minus thin curve
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no {series} marker found near x={xv:g}"})
            continue
        # Keep the largest connected component inside the bin slice
        lab, n = ndimage.label(sub)
        sizes = [int((lab == i).sum()) for i in range(1, n + 1)]
        best_i = int(np.argmax(sizes)) + 1
        br, bc = np.nonzero(lab == best_i)
        r_min, r_max = float(br.min()) + 20.0, float(br.max()) + 20.0
        c_min, c_max = float(bc.min()) + c0, float(bc.max()) + c0
        # For squares/circles: bounding-box midpoint == centroid.
        # For OriginLab up/down triangles: center of bounding box (0.5*(r_min+r_max)) aligns with plotted data value!
        # For partially occluded gray square (series 005) where red covers the bottom or top: check height vs expected ~65 px!
        h_blob = r_max - r_min + 1.0
        if series == "005" and h_blob < 52.0:
            # Partially occluded square: use visible centroid with wider uncertainty covering occlusion
            r_mid = float(br.mean()) + 20.5
            err_px = max(6.0, 0.5 * (65.0 - h_blob)) + y_cal.res_px
        else:
            r_mid = 0.5 * (r_min + r_max) + 0.5
            err_px = 4.0 + y_cal.res_px
        c_mid = 0.5 * (c_min + c_max) + 0.5
        y_val = y_cal.to_val(r_mid)
        ya = y_cal.to_val(r_mid - err_px)
        yb = y_cal.to_val(r_mid + err_px)
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": min(ya, yb),
            "hi": max(ya, yb),
            "px": [c_mid, r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"marker {series} ({len(br)} px, h={h_blob:.0f}px), rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
        })
    return out


def read_fig2d_fit(image: np.ndarray, series: str, at: list[float]) -> list[dict]:
    """Read smooth log-normal fit curve in Fig. 2D (stripping solid markers so only the curve stroke is read)."""
    im = image.astype(float)
    H, W, _ = im.shape
    hr, vl, x_cal, y_cal = _calibrate_fig2d(image)

    mask = _series_mask(im, series)
    mask[: int(y_cal.to_px(0.37)), int(x_cal.to_px(2.75)) :] = False
    # Remove thick solid markers (opening 22x22 dilated slightly) so curve readings aren't pulled by markers!
    markers = ndimage.binary_opening(mask, structure=np.ones((22, 22), dtype=bool))
    markers_dil = ndimage.binary_dilation(markers, structure=np.ones((9, 9), dtype=bool))
    curve_only = mask & (~markers_dil)

    # Note: where the curve passes right behind a marker of its own color, `at` points chosen between 0.2 um bins
    # (e.g. 0.3, 0.5, 0.7, 0.9, ...) hit the pure unobstructed curve!
    r_top = 20
    r_bot = hr[0] - 2
    c_left = vl[1] + 5
    c_right = W - 10
    return read_colored_curve_at(
        curve_only, x_cal, y_cal, at, r_top, r_bot, c_left, c_right,
        search_cols=12, window_half=3,
    )
