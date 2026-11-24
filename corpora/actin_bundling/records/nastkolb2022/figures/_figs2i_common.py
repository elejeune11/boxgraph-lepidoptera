"""Shared calibration and 1 uM VASP histogram/log-normal reader for Fig. S2I (PDF p. 16, 1200 dpi)."""

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


def _calibrate_figs2i(image: np.ndarray):
    im = image.astype(float)
    H, W, _ = im.shape
    gray = im.mean(axis=2)
    dark = gray < 105.0
    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: 0, 2, 4, 6, 8 um (5 major ticks)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=48, min_len=12, coord_min=vl[1], coord_max=W)
    xt = select_major_ticks(xt_raw, expected=5, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 2.0, 4.0, 6.0, 8.0], max_res_px=2.5)

    # Y axis: 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0 (7 major ticks)
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=48, min_len=12, coord_min=50, coord_max=hr[0] - 5)
    yt = select_major_ticks(yt_raw, expected=7, frac=0.72)
    y_cal = LinearCal(yt, [0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0], max_res_px=2.5)
    return hr, vl, x_cal, y_cal


def _purple_mask(im: np.ndarray) -> np.ndarray:
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    return (R > 140) & (B > 175) & (G < 155) & (B - G > 50) & (R - G > 25)


def read_figs2i_hist_100(image: np.ndarray, bin_centers: list[float]) -> list[dict]:
    """Read 1 uM VASP purple diamond histogram markers (0.2..8.4 um, 0.2 um bins) in Fig. S2I."""
    im = image.astype(float)
    H, W, _ = im.shape
    hr, vl, x_cal, y_cal = _calibrate_figs2i(image)

    mask = _purple_mask(im)
    # Exclude legend box (x > 5.0, y > 0.28)
    mask[: int(y_cal.to_px(0.27)), int(x_cal.to_px(4.9)) :] = False
    mask[hr[0] - 2 :, :] = False
    mask[:, : vl[1] + 4] = False

    # Note: 0.2 um = 52.3 px, diamonds are ~55 px wide so adjacent diamonds touch horizontally!
    # Instead of opening + connected components, slice each exact 0.2 um bin column +/- 15 px around bin center!
    markers = ndimage.binary_opening(mask, structure=np.ones((14, 14), dtype=bool))
    out = []
    half_w = 16
    for xv in bin_centers:
        pc = x_cal.to_px(xv)
        c0 = max(vl[1] + 4, int(round(pc)) - half_w)
        c1 = min(W - 4, int(round(pc)) + half_w)
        sub = markers[int(y_cal.to_px(0.20)) : hr[0] - 2, c0:c1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 40:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no purple diamond near x={xv:g}"})
            continue
        r0 = int(y_cal.to_px(0.20))
        # Take exact vertical bounding-box midpoint at the central 10 columns around pc (diamond vertices align at pc!)
        narrow = markers[r0 : hr[0] - 2, max(c0, int(round(pc)) - 5) : min(c1, int(round(pc)) + 6)]
        nrr, _ = np.nonzero(narrow)
        if len(nrr) >= 20:
            r_mid = 0.5 * (float(nrr.min()) + float(nrr.max())) + 0.5 + r0
        else:
            r_mid = 0.5 * (float(rr.min()) + float(rr.max())) + 0.5 + r0
        y_val = y_cal.to_val(r_mid)
        err_px = 5.0 + y_cal.res_px
        ya, yb = y_cal.to_val(r_mid - err_px), y_cal.to_val(r_mid + err_px)
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": min(ya, yb),
            "hi": max(ya, yb),
            "px": [pc, r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"1uM purple diamond at x={xv:g}, rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
        })
    return out


def read_figs2i_fit_100(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read 1 uM VASP purple log-normal fit curve in Fig. S2I."""
    im = image.astype(float)
    H, W, _ = im.shape
    hr, vl, x_cal, y_cal = _calibrate_figs2i(image)

    mask = _purple_mask(im)
    mask[: int(y_cal.to_px(0.27)), int(x_cal.to_px(4.9)) :] = False
    # Detect each 0.2 um diamond center and erase its exact L1 diamond footprint (|dr|+|dc| <= 29 px)
    # so all 4 pointy diamond tips are cleanly removed while preserving the smooth fit curve!
    bins_all = [round(0.2 * i, 1) for i in range(1, 43)]
    hist_pts = read_figs2i_hist_100(image, bins_all)
    curve_only = mask.copy()
    R_dia = 30
    yy, xx = np.ogrid[-R_dia : R_dia + 1, -R_dia : R_dia + 1]
    dia_mask = (np.abs(yy) + np.abs(xx)) <= R_dia
    for pt in hist_pts:
        if pt.get("y") is None:
            continue
        ci = int(round(pt["px"][0]))
        ri = int(round(pt["px"][1]))
        r0, r1 = max(0, ri - R_dia), min(H, ri + R_dia + 1)
        c0, c1 = max(0, ci - R_dia), min(W, ci + R_dia + 1)
        curve_only[r0:r1, c0:c1] &= ~dia_mask[r0 - (ri - R_dia) : r1 - (ri - R_dia), c0 - (ci - R_dia) : c1 - (ci - R_dia)]
    return read_colored_curve_at(
        curve_only, x_cal, y_cal, at, int(y_cal.to_px(0.18)), hr[0] - 2, vl[1] + 5, W - 5,
        search_cols=6, window_half=2,
    )
