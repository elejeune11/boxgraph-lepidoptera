"""Shared calibrated reader for Fig. 1C (FRAP normalized fluorescence intensity vs time)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    xt = C.match_ticks_regular(xt, 4, tol=4.0)
    yt = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(8, 45))
    yt = C.match_ticks_regular(yt, 4, tol=4.0)

    xa = C.verify_axis("linear", list(zip(xt, [0.0, 5.0, 10.0, 15.0])) + [(vr0[0], 20.0)], max_residual_px=1.5)
    ya = C.verify_axis("linear", [(ht[0], 1.0)] + list(zip(yt, [0.8, 0.6, 0.4, 0.2])) + [(hb[0], 0.0)], max_residual_px=1.5)
    return dark, xa, ya, int(ht[2]) + 68, int(hb[1]) - 68, int(vl0[2]) + 68, int(vr0[1]) - 68


def _marker_mask(dark: np.ndarray, r0: int, r1: int, c0: int, c1: int):
    m = np.zeros_like(dark)
    m[r0:r1, c0:c1] = dark[r0:r1, c0:c1]
    # Fast distance transform separates solid circles (radius ~44 px) from thin dashed line (half-width ~7 px)
    dist = ndimage.distance_transform_edt(m)
    cores = dist >= 22.0
    return m, cores


def read_data(image: np.ndarray) -> list[dict]:
    dark, xa, ya, r0, r1, c0, c1 = calibrate(image)
    _, cores = _marker_mask(dark, r0, r1, c0, c1)
    lab, n = ndimage.label(cores)
    out = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 200:
            continue
        c_m = float(cc.mean()) + 0.5
        r_m = float(rr.mean()) + 0.5
        half_px = 3.0 + ya.residual_px
        yv = float(ya.to_value(r_m))
        ya_, yb_ = float(ya.to_value(r_m - half_px)), float(ya.to_value(r_m + half_px))
        out.append(
            {
                "x": float(xa.to_value(c_m)),
                "y": yv,
                "lo": min(ya_, yb_),
                "hi": max(ya_, yb_),
                "px": [c_m, r_m],
                "px_lo": float(r_m - half_px),
                "px_hi": float(r_m + half_px),
                "note": f"Fig1C circle core={len(rr)}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            }
        )
    out.sort(key=lambda d: d["x"])
    return out


def read_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, xa, ya, r0, r1, c0, c1 = calibrate(image)
    m, cores = _marker_mask(dark, r0, r1, c0, c1)
    # Erase every circular marker with generous margin (dot radius ~46 px -> mask out radius 66 px)
    dash = m.copy()
    H, W = dash.shape
    for b in C.find_blobs(cores, min_area=200):
        cm, rm = b["c"] - 0.5, b["r"] - 0.5
        r_lo, r_hi = max(0, int(rm) - 70), min(H, int(rm) + 71)
        c_lo, c_hi = max(0, int(cm) - 70), min(W, int(cm) + 71)
        yy, xx = np.ogrid[r_lo:r_hi, c_lo:c_hi]
        dash[r_lo:r_hi, c_lo:c_hi] &= ((xx - cm) ** 2 + (yy - rm) ** 2) > 66.0**2

    # Keep only genuine thin dash strokes (drop any tiny stray specks < 80 px)
    lab, _ = ndimage.label(dash)
    clean = np.zeros_like(dash)
    for b in C.find_blobs(dash, min_area=80, max_area=4000):
        if b["h"] <= 70 and b["w"] <= 140:
            sl = (slice(b["r0"], b["r1"] + 1), slice(b["c0"], b["c1"] + 1))
            clean[sl] |= lab[sl] == b["label"]

    return [
        C.sample_dashed_curve_y(
            clean, xa, ya, float(x), r0, r1, max_search_dx=320, sigma_dx=140.0, deg=2, note_prefix="Fig1C exp fit: "
        )
        for x in at
    ]
