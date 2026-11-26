"""Shared calibration and marker detection for Fig. 2b (p. 2 / printed p. 790)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy import ndimage
import _common as C


def read_series(image: np.ndarray, series: str) -> list[dict]:
    box = C.find_plot_box(image, dark_thr=120.0, sat_thr=35.0, min_h_frac=0.6, min_v_frac=0.5)
    dark = C.dark_neutral_mask(image, dark_thr=120.0, sat_thr=35.0)
    top, bot, left, right = box["top"], box["bot"], box["left"], box["right"]

    xt_all = C.scan_ticks(dark, "x", bot[0] - 26, bot[0] - 6, left[1] + 5, right[0] - 5, min_frac=0.6, max_len=20)
    yt_all = C.scan_ticks(dark, "y", top[1] + 5, bot[0] - 5, left[1] + 6, left[1] + 26, min_frac=0.6, max_len=20)

    xt_reg = C.select_regular_ticks(xt_all, 8, tol=2.5)
    xa = C.make_axis("linear", list(zip(xt_reg, [280.0, 330.0, 380.0, 430.0, 480.0, 530.0, 580.0, 630.0])), max_res_px=1.5)

    yt_reg = C.select_regular_ticks(yt_all, 7, tol=2.5)
    ya = C.make_axis("linear", list(zip(yt_reg, [0.20, 0.15, 0.10, 0.05, 0.00, -0.05, -0.10])), max_res_px=1.5)

    rgb = image.astype(np.float32)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    interior = np.zeros(image.shape[:2], dtype=bool)
    interior[top[1] + 8 : bot[0] - 8, left[1] + 8 : right[0] - 8] = True

    blue_m = interior & (B > 110.0) & ((B - R) > 45.0) & ((B - G) > 35.0)
    red_m = interior & (R > 140.0) & ((R - G) > 60.0) & ((R - B) > 60.0)
    black_m = interior & (rgb.max(axis=2) < 90.0)

    def detect_circles(m: np.ndarray) -> list[tuple[float, float, int]]:
        # Fast O(N) Euclidean distance transform separates circular markers (radius ~38 px) from connecting lines (half-width <= 18 px)
        core = ndimage.distance_transform_edt(m) >= 19.5
        lab, n = ndimage.label(core)
        pts = []
        for i in range(1, n + 1):
            rr, cc = np.nonzero(lab == i)
            if len(rr) < 120:
                continue
            pts.append((float(cc.mean()) + 0.5, float(rr.mean()) + 0.5, int(len(rr))))
        pts.sort(key=lambda t: t[0])
        return pts

    blue_pts = detect_circles(blue_m)
    if len(blue_pts) != 5:
        raise RuntimeError(f"Fig. 2b: expected 5 blue tip circles, found {len(blue_pts)}")

    if series == "tip":
        pts = [(c, r, f"Fig. 2b tip (blue circle, area={a}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px)", 2.0)
               for c, r, a in blue_pts]
    elif series in ("proximal", "distal"):
        m = black_m if series == "proximal" else red_m
        c0, r0, _ = blue_pts[0]
        det = [p for p in detect_circles(m) if p[0] > c0 + 60.0]
        if len(det) != 4:
            raise RuntimeError(f"Fig. 2b {series}: expected 4 unoccluded circles after t0, found {len(det)}")
        # All three trajectories originate at the shared t0 ~ 280 s, Delta L = 0 reference marker
        # (blue tip circle is drawn on top of the red and black circles at t0).
        label = "proximal bleach edge (black)" if series == "proximal" else "distal bleach edge (red)"
        pts = [(c0, r0, f"Fig. 2b {label}: shared t0 origin (occluded by blue marker), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px", 6.0)]
        for c, r, a in det:
            pts.append((c, r, f"Fig. 2b {label} (circle, area={a}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px)", 2.5))
    else:
        raise ValueError(f"Unknown series {series!r}")

    out = []
    for c, r, note, half_px in pts:
        out.append(C.reading_point(x_val=xa.to_value(c), col=c, row=r, y_axis=ya, half_px=half_px, note=note))
    return out
