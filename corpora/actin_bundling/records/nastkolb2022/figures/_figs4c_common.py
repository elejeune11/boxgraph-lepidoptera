"""Shared calibration, histogram marker reader, and log-normal fit curve reader for Fig. S4C (PDF p. 18, 1200 dpi)."""

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
)


def _calibrate_figs4c(image: np.ndarray):
    im = image.astype(float)
    H, W, _ = im.shape
    gray = im.mean(axis=2)
    dark = gray < 105.0
    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 um (7 major ticks)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=55, min_len=14, coord_min=vl[1] + 50, coord_max=W)
    xt = select_major_ticks(xt_raw, expected=7, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], max_res_px=2.5)

    # Y axis: 0.5, 0.4, 0.3, 0.2, 0.1, 0.0 (6 major ticks)
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=55, min_len=14, coord_min=40, coord_max=hr[0] - 5)
    yt = select_major_ticks(yt_raw, expected=6, frac=0.72)
    y_cal = LinearCal(yt, [0.5, 0.4, 0.3, 0.2, 0.1, 0.0], max_res_px=2.5)
    return hr, vl, x_cal, y_cal


def _s4c_mask(im: np.ndarray, series: str) -> np.ndarray:
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    if series == "vasp":
        # Blue up-triangles (~RGB 20,105,220)
        return (B > 155) & (B - R > 100) & (B - G > 40)
    elif series == "nocp_slb":
        # Olive down-triangles (~RGB 128,128,0)
        return (R > 95) & (R < 175) & (G > 95) & (G < 175) & (B < 45) & (np.abs(R - G) < 35)
    elif series == "nocp_ml":
        # Orange circles (~RGB 245,100,0)
        return (R > 195) & (G > 60) & (G < 150) & (B < 45) & (R - G > 85)
    else:
        raise ValueError(series)


def read_figs4c_hist(image: np.ndarray, series: str, bin_centers: list[float]) -> list[dict]:
    """Read discrete histogram markers at exact 0.2 um bin centers in Fig. S4C."""
    im = image.astype(float)
    H, W, _ = im.shape
    hr, vl, x_cal, y_cal = _calibrate_figs4c(image)

    mask = _s4c_mask(im, series)
    # Exclude legend box (top right: x > 1.4, y > 0.31)
    mask[: int(y_cal.to_px(0.30)), int(x_cal.to_px(1.38)) :] = False
    mask[hr[0] - 2 :, :] = False
    mask[:, : vl[1] + 5] = False

    markers = ndimage.binary_opening(mask, structure=np.ones((12, 12), dtype=bool))
    out = []
    half_bin_px = int(0.08 * abs(x_cal.to_px(1.0) - x_cal.to_px(0.0)))
    for xv in bin_centers:
        pc = x_cal.to_px(xv)
        c0 = max(vl[1] + 5, int(round(pc)) - half_bin_px)
        c1 = min(W - 5, int(round(pc)) + half_bin_px)
        sub = markers[40 : hr[0] - 2, c0:c1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 40:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no {series} marker near x={xv:g}"})
            continue
        lab, n = ndimage.label(sub)
        sizes = [int((lab == i).sum()) for i in range(1, n + 1)]
        best_i = int(np.argmax(sizes)) + 1
        br, bc = np.nonzero(lab == best_i)
        r_mid = 0.5 * (float(br.min()) + float(br.max())) + 40.5
        c_mid = 0.5 * (float(bc.min()) + float(bc.max())) + c0 + 0.5
        y_val = y_cal.to_val(r_mid)
        err_px = 4.5 + y_cal.res_px
        ya, yb = y_cal.to_val(r_mid - err_px), y_cal.to_val(r_mid + err_px)
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": min(ya, yb),
            "hi": max(ya, yb),
            "px": [c_mid, r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"S4C marker {series} ({len(br)} px), rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
        })
    return out


def read_figs4c_fit(image: np.ndarray, series: str, at: list[float]) -> list[dict]:
    """Read smooth log-normal fit curves in Fig. S4C (stripping solid markers)."""
    im = image.astype(float)
    H, W, _ = im.shape
    hr, vl, x_cal, y_cal = _calibrate_figs4c(image)

    mask = _s4c_mask(im, series)
    mask[: int(y_cal.to_px(0.30)), int(x_cal.to_px(1.38)) :] = False
    markers = ndimage.binary_opening(mask, structure=np.ones((20, 20), dtype=bool))
    curve_only = mask & (~ndimage.binary_dilation(markers, structure=np.ones((7, 7), dtype=bool)))
    return read_colored_curve_at(
        curve_only, x_cal, y_cal, at, 40, hr[0] - 2, vl[1] + 5, W - 5,
        search_cols=10, window_half=3,
    )
