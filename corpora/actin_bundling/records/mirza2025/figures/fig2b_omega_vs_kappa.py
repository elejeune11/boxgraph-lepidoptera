"""Calibrated reader for Fig. 2(b) (PDF p. 36): architectural order parameter omega vs kappa.
Detects plot box, 9 inward x-ticks (-0.8..0.8), 4 inward y-ticks (-0.04, -0.02, 0.02, 0.04),
and all 18 dark square simulation markers directly from the rendered 1200 dpi crop.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    frame = C.find_box(image, dark=160, max_sat=80, min_frac=0.55)
    xt_c = C.detect_x_ticks(image, frame, side="in", spine="bottom", dark=160, max_sat=80, min_len=30, max_len=150)
    yt_c = C.detect_y_ticks(image, frame, side="in", spine="left", dark=160, max_sat=80, min_len=30, max_len=150)

    xt = C.select_regular_ticks(xt_c, expected_n=9, tol_px=12.0)
    xa = C.calibrate_axis(xt, [-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8], max_residual_px=5.0)

    # Filter inward y-ticks by stroke width (excluding left-border marker touching spine at kappa=-0.8)
    yt_inner = sorted(p for p, l, w in yt_c if 40 <= l <= 80 and w <= 8)
    ya = C.calibrate_axis(yt_inner, [0.04, 0.02, -0.02, -0.04], max_residual_px=4.0)

    # Detect solid dark square markers (~48x48 px at 1200 dpi) via fast separable box filter
    dn = C.neutral_dark_mask(image, dark=120, max_sat=80)
    c0, r0, c1, r1 = [int(v) for v in frame["box"]]
    inner = np.zeros_like(dn, dtype=np.float32)
    inner[r0 + 5 : r1 - 5, max(0, c0 - 25) : min(dn.shape[1], c1 + 25)] = dn[
        r0 + 5 : r1 - 5, max(0, c0 - 25) : min(dn.shape[1], c1 + 25)
    ]
    cores = ndimage.uniform_filter(inner, size=27) > 0.85
    lab, n = ndimage.label(cores)

    # Native raster resolution is ~300 dpi (4 px per native pixel at 1200 dpi)
    native_px = 4.0
    h_px = 0.5 * native_px + max(ya.residual_px, 0.5 * native_px) + 1.5

    rows = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 100:
            continue
        r_m = float(rr.mean() + 0.5)
        c_m = float(cc.mean() + 0.5)
        xv = xa.to_value(c_m)
        yv = ya.to_value(r_m)
        a, b = ya.to_value(r_m - h_px), ya.to_value(r_m + h_px)
        rows.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [c_m, r_m],
            "px_lo": r_m - h_px,
            "px_hi": r_m + h_px,
            "note": f"square marker core={len(rr)}px cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
        })
    rows.sort(key=lambda d: d["x"])
    return rows
