"""Shared calibration and 6-series kinetics reader for Fig. 3E (PDF p. 5, 1200 dpi)."""

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


def read_fig3e(image: np.ndarray, series: str, at: list[float]) -> list[dict]:
    """Read Fig. 3E TIRFM actin polymerization kinetics (6 distinct colored marker series)."""
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 min (11 major ticks)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=42, min_len=14, coord_min=vl[1] + 50, coord_max=W - 10)
    xt = select_major_ticks(xt_raw, expected=11, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0], max_res_px=3.0)

    # Y axis: 1.0, 0.8, 0.6, 0.4, 0.2, 0.0 (6 major ticks)
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=45, min_len=14, coord_min=100, coord_max=hr[0] - 20)
    yt = select_major_ticks(yt_raw, expected=6, frac=0.65)
    y_cal = LinearCal(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0], max_res_px=3.0)

    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    if series == "vasp":
        # Blue up-triangles (~RGB 20,105,220)
        mask = (B > 155) & (B - R > 100) & (B - G > 40)
    elif series == "vasp_dtetra":
        # Gold/olive right-triangles (~RGB 200,148,0)
        mask = (R > 155) & (G > 110) & (G < 185) & (B < 55) & (R - G > 20) & (R - G < 85)
    elif series == "actinin_050":
        # Light gray diamonds (~RGB 192,192,192) along bottom y < 0.16
        mask = (gray > 160) & (gray < 215) & (sat < 16)
        mask[: int(y_cal.to_px(0.18)), :] = False
    elif series == "actinin_005_nocp":
        # Cyan/teal left-triangles (~RGB 0,195,195)
        mask = (G > 150) & (B > 150) & (R < 55) & (np.abs(G - B) < 40)
    elif series == "mdia1":
        # Dark brown hexagons (~RGB 120,72,72)
        mask = (R > 85) & (R < 155) & (G > 45) & (G < 105) & (B > 45) & (B < 105) & (R - G > 25) & (np.abs(G - B) < 25)
    elif series == "mdia1_actinin":
        # Bright orange pentagons (~RGB 245,95,0)
        mask = (R > 200) & (G > 55) & (G < 140) & (B < 45) & (R - G > 90)
    else:
        raise ValueError(series)

    # Mask out legend box: x in [9.5, 20.5], y in [0.16, 0.62]
    mask[int(y_cal.to_px(0.62)) : int(y_cal.to_px(0.16)), int(x_cal.to_px(9.4)) :] = False
    mask[hr[0] - 2 :, :] = False
    mask[:, : vl[1] + 4] = False

    # Clean tiny anti-aliasing specks
    clean = ndimage.binary_opening(mask, structure=np.ones((5, 5), dtype=bool))
    px_per_min = abs(x_cal.to_px(2.0) - x_cal.to_px(0.0)) / 2.0

    # If `at` is empty: detect every discrete connected-component marker (great for sparse mDia1!)
    # If `at` is provided: at each requested time t, find the centroid of clean pixels within +/- 0.45 min (or nearest blob for sparse series!)
    out = []
    if not at:
        lab, n = ndimage.label(clean)
        pts = []
        for k in range(1, n + 1):
            rr, cc = np.nonzero(lab == k)
            if len(rr) < 120:
                continue
            r_mid = 0.5 * (float(rr.min()) + float(rr.max())) + 0.5
            c_mid = 0.5 * (float(cc.min()) + float(cc.max())) + 0.5
            pts.append((c_mid, r_mid, len(rr)))
        pts.sort(key=lambda p: p[0])
        for c_mid, r_mid, area in pts:
            xv = x_cal.to_val(c_mid)
            yv = y_cal.to_val(r_mid)
            err_px = 4.0 + y_cal.res_px
            ya, yb = y_cal.to_val(r_mid - err_px), y_cal.to_val(r_mid + err_px)
            out.append({
                "x": round(xv, 2),
                "y": yv,
                "lo": min(ya, yb),
                "hi": max(ya, yb),
                "px": [c_mid, r_mid],
                "px_lo": r_mid - err_px,
                "px_hi": r_mid + err_px,
                "note": f"marker {series} ({area} px), rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
            })
        return out

    # Sampled at requested times `at`: find nearest marker pixels within +/- 0.85 min (narrower at early steep rise)
    for tv in at:
        pc = x_cal.to_px(tv)
        half_w = int((0.25 if tv < 5.0 else 0.75) * px_per_min)
        c0 = max(vl[1] + 5, int(round(pc)) - half_w)
        c1 = min(W - 5, int(round(pc)) + half_w)
        sub = clean[40 : hr[0] - 2, c0:c1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 50:
            out.append({"x": float(tv), "y": None, "found": False, "note": f"no {series} marker near t={tv:g}"})
            continue
        if tv < 5.0:
            # Dense rising region: evaluate via local linear fit (or narrow band +/- 10 px around pc)
            narrow = clean[40 : hr[0] - 2, max(vl[1] + 5, int(round(pc)) - 12) : min(W - 5, int(round(pc)) + 13)]
            nrr, ncc = np.nonzero(narrow)
            if len(nrr) >= 20:
                r_mid = float(nrr.mean()) + 40.5
                c_mid = pc
                t_actual = float(tv)
            else:
                r_mid = float(rr.mean()) + 40.5
                c_mid = float(cc.mean()) + c0 + 0.5
                t_actual = round(x_cal.to_val(c_mid), 2)
            err_px = 8.0 + y_cal.res_px
        else:
            # Discrete marker region: pick connected component closest in x to pc!
            lab, n = ndimage.label(sub)
            best = None
            for k in range(1, n + 1):
                kr, kc = np.nonzero(lab == k)
                if len(kr) < 80:
                    continue
                cm = 0.5 * (float(kc.min()) + float(kc.max())) + 0.5 + c0
                rm = 0.5 * (float(kr.min()) + float(kr.max())) + 40.5
                if best is None or abs(cm - pc) < abs(best[0] - pc):
                    best = (cm, rm, len(kr))
            if best is None:
                out.append({"x": float(tv), "y": None, "found": False, "note": f"no {series} blob near t={tv:g}"})
                continue
            c_mid, r_mid, _ = best
            t_actual = round(x_cal.to_val(c_mid), 2)
            err_px = 5.0 + y_cal.res_px
        yv = y_cal.to_val(r_mid)
        ya, yb = y_cal.to_val(r_mid - err_px), y_cal.to_val(r_mid + err_px)
        out.append({
            "x": t_actual,
            "y": yv,
            "lo": min(ya, yb),
            "hi": max(ya, yb),
            "px": [c_mid, r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"{series} at t={t_actual:g}min (req {tv:g}), rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
        })
    return out
