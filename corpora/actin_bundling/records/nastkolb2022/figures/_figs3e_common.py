"""Shared calibration and 4-series kinetics reader for Fig. S3E actin titration (PDF p. 17, 1200 dpi)."""

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


def read_figs3e(image: np.ndarray, series: str, at: list[float]) -> list[dict]:
    """Read Fig. S3E raw TIRFM actin intensity kinetics at 0.5, 1, 2, 4 uM actin."""
    im = image.astype(float)
    H, W, _ = im.shape
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    gray = im.mean(axis=2)
    dark = gray < 105.0

    h_spines, v_spines = find_black_spines(gray, dark=105.0, min_frac_h=0.45, min_frac_v=0.45)
    hr = h_spines[-1]
    vl = v_spines[0]

    # X axis: 0, 10, 20, 30, 40, 50, 60 min (7 major ticks)
    xt_raw = detect_ticks_along_spine(dark, hr, "x", "down", min_offset=2, max_offset=48, min_len=14, coord_min=vl[1] + 20, coord_max=W)
    xt = select_major_ticks(xt_raw, expected=7, frac=0.72)
    x_cal = LinearCal(xt, [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0], max_res_px=2.5)

    # Y axis: 6000, 5000, 4000, 3000, 2000, 1000, 0 AU (7 major ticks)
    yt_raw = detect_ticks_along_spine(dark, vl, "y", "left", min_offset=2, max_offset=48, min_len=12, coord_min=40, coord_max=hr[0] - 5)
    yt = select_major_ticks(yt_raw, expected=7, frac=0.72)
    y_cal = LinearCal(yt, [6000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0, 0.0], max_res_px=2.5)

    sat = np.maximum(np.maximum(R, G), B) - np.minimum(np.minimum(R, G), B)
    if series == "05":
        mask = (gray > 50) & (gray < 125) & (sat < 22)
    elif series == "1":
        mask = (R > 180) & (R - G > 110) & (R - B > 110)
    elif series == "2":
        mask = (B > 155) & (B - R > 110) & (B - G > 40)
    elif series == "4":
        mask = (G > 125) & (G - R > 60) & (G - B > 25)
    else:
        raise ValueError(series)

    # Mask out legend box (top left: x < 18, y > 4450)
    mask[: int(y_cal.to_px(4450.0)), : int(x_cal.to_px(18.0))] = False
    mask[hr[0] - 3 :, :] = False
    mask[:, : vl[1] + 5] = False

    clean = ndimage.binary_opening(mask, structure=np.ones((6, 6), dtype=bool))
    px_per_min = abs(x_cal.to_px(10.0) - x_cal.to_px(0.0)) / 10.0

    out = []
    for tv in at:
        pc = x_cal.to_px(tv)
        half_w = int((0.35 if tv < 7.0 else 1.4) * px_per_min)
        c0 = max(vl[1] + 5, int(round(pc)) - half_w)
        c1 = min(W - 5, int(round(pc)) + half_w)
        sub = clean[40 : hr[0] - 3, c0:c1]
        rr, cc = np.nonzero(sub)
        if len(rr) < 40:
            out.append({"x": float(tv), "y": None, "found": False, "note": f"no {series} marker near t={tv:g}"})
            continue
        if tv < 7.0:
            narrow = clean[40 : hr[0] - 3, max(vl[1] + 5, int(round(pc)) - 8) : min(W - 5, int(round(pc)) + 9)]
            nrr, _ = np.nonzero(narrow)
            if len(nrr) >= 20:
                r_mid = float(nrr.mean()) + 40.5
                c_mid = pc
                t_actual = float(tv)
            else:
                r_mid = float(rr.mean()) + 40.5
                c_mid = float(cc.mean()) + c0 + 0.5
                t_actual = round(x_cal.to_val(c_mid), 2)
            err_px = 10.0 + y_cal.res_px
        else:
            lab, n = ndimage.label(sub)
            best = None
            for k in range(1, n + 1):
                kr, kc = np.nonzero(lab == k)
                if len(kr) < 80:
                    continue
                cm = 0.5 * (float(kc.min()) + float(kc.max())) + 0.5 + c0
                rm = 0.5 * (float(kr.min()) + float(kr.max())) + 40.5
                if best is None or abs(cm - pc) < abs(best[0] - pc):
                    best = (cm, rm)
            if best is None:
                out.append({"x": float(tv), "y": None, "found": False, "note": f"no {series} blob near t={tv:g}"})
                continue
            c_mid, r_mid = best
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
            "note": f"actin {series}uM at t={t_actual:g}min (req {tv:g}), rx={x_cal.res_px:.2f}px, ry={y_cal.res_px:.2f}px",
        })
    return out
