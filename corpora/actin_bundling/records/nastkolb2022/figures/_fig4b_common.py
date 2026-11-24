"""Shared calibration and 4-series marker reader for Fig. 4B Spearman colocalization vs time (PDF p. 7, 1200 dpi)."""

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
    select_major_ticks,
)


def read_fig4b(image: np.ndarray, series: str) -> list[dict]:
    """Read all discrete plotted markers for one of the 4 series in Fig. 4B."""
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: 0, 10, 20, 30, 40 min (5 major ticks)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=55, min_len=14, coord_min=vl[1] + 50, coord_max=W - 10)
    xt = select_major_ticks(xt_raw, expected=5, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 10.0, 20.0, 30.0, 40.0], max_res_px=2.0)

    # Y axis: 1.0, 0.8, 0.6, 0.4, 0.2, 0.0 (6 major ticks)
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=55, min_len=14, coord_min=20, coord_max=hr[1] + 10)
    yt = select_major_ticks(yt_raw, expected=6, frac=0.72)
    y_cal = LinearCal(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0], max_res_px=3.0)

    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    if series == "vasp_slb":
        # Black squares
        mask = (R < 55) & (G < 55) & (B < 55)
    elif series == "vasp_ml":
        # Red up-triangles
        mask = (R > 185) & (R - G > 110) & (R - B > 110)
    elif series == "dtetra_slb":
        # Blue diamonds
        mask = (B > 155) & (B - R > 110) & (B - G > 40)
    elif series == "dtetra_ml":
        # Green left-triangles
        mask = (G > 125) & (G - R > 60) & (G - B > 25)
    else:
        raise ValueError(series)

    # Mask out legend box (x in [25, 44], y in [0.42, 0.67]) and axes borders
    mask[int(y_cal.to_px(0.68)) : int(y_cal.to_px(0.41)), int(x_cal.to_px(24.5)) : int(x_cal.to_px(44.5))] = False
    mask[hr[0] - 5 :, :] = False
    mask[:, : vl[1] + 10] = False
    mask[:30, :] = False

    clean = ndimage.binary_opening(mask, structure=np.ones((12, 12), dtype=bool))
    lab, n = ndimage.label(clean)

    pts = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 400:
            continue
        # OriginLab symbols are centered at the bounding-box midpoint (note: at t=0 and t=11 green left-triangle
        # slightly overlaps red up-triangle horizontally; vertical bounding-box midpoint is still accurate!)
        r_mid = 0.5 * (float(rr.min()) + float(rr.max())) + 0.5
        c_mid = 0.5 * (float(cc.min()) + float(cc.max())) + 0.5
        pts.append((c_mid, r_mid, len(rr)))
    pts.sort(key=lambda p: p[0])

    out = []
    for c_mid, r_mid, area in pts:
        xv = x_cal.to_val(c_mid)
        yv = y_cal.to_val(r_mid)
        err_px = 4.0 + y_cal.res_px
        ya, yb = y_cal.to_val(r_mid - err_px), y_cal.to_val(r_mid + err_px)
        out.append({
            "x": round(max(0.0, xv), 2) if abs(xv) < 0.3 else round(xv, 2),
            "y": yv,
            "lo": min(ya, yb),
            "hi": max(ya, yb),
            "px": [c_mid, r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"marker {series} ({area} px), rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
        })
    return out
