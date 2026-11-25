"""Calibrated reader for Fig. 4B (top panel, left y-axis, black; PDF p. 6 / printed p. 26575):
Lo lipid domain area fraction phi_Lo vs. time t [min] during cholesterol insertion (t = 0 to 11 min).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from boxgraph import digitize as D


def _group(indices: np.ndarray, gap: int = 2) -> list[list[int]]:
    out: list[list[int]] = []
    for x in indices:
        xi = int(x)
        if out and xi <= out[-1][-1] + gap:
            out[-1].append(xi)
        else:
            out.append([xi])
    return out


def _calibrate(image: np.ndarray):
    H, W = image.shape[:2]
    r, g, b = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    black = (mx < 130) & ((mx - mn) < 35)

    # Two horizontal black x-axes (top panel and bottom panel)
    hg = _group(np.where(black.sum(axis=1) > 0.5 * W)[0], gap=6)
    r_top_x0, r_top_x1 = hg[0][0], hg[0][-1]
    r_bot_x0, r_bot_x1 = hg[1][0], hg[1][-1]

    # Top panel left black spine (search columns 600..1000 to skip y-axis labels)
    tl_g = _group(np.where(black[:r_top_x0, 600:1000].sum(axis=0) > 0.5 * r_top_x0)[0] + 600, gap=4)[0]
    c_tl0, c_tl1 = tl_g[0], tl_g[-1]
    c_zero = (c_tl0 + c_tl1) / 2.0 + 0.5

    # Bottom x-axis ticks (5, 10, 15, 20 min) share exact columns with top panel
    xt_raw = []
    for c in range(c_tl1 + 20, W - 10):
        n = 0
        while r_bot_x0 - 1 - n >= 0 and black[r_bot_x0 - 1 - n, c]:
            n += 1
        if 14 <= n <= 60:
            xt_raw.append(c)
    xt_groups = [g for g in _group(np.array(xt_raw), gap=1) if 8 <= len(g) <= 30]
    xt_pos = sorted((g[0] + g[-1]) / 2.0 + 0.5 for g in xt_groups)
    x_ticks = D.regular_ticks([c_zero] + xt_pos, tol=4.0)
    if len(x_ticks) != 5:
        raise D.CalibrationError(f"Fig. 4B top: expected 5 x-ticks (0,5,10,15,20), got {x_ticks}")
    xa = D.Axis("linear", list(zip(x_ticks, [0.0, 5.0, 10.0, 15.0, 20.0])))

    # Top panel left y-axis ticks (0.50, 0.45, 0.40, 0.35, 0.30; pointing rightward from left spine)
    yt_raw = []
    for r_idx in range(100, r_top_x0 - 10):
        n = 0
        while c_tl1 + 1 + n < W and black[r_idx, c_tl1 + 1 + n]:
            n += 1
        if 25 <= n <= 80:
            yt_raw.append(r_idx)
    yt_groups = [g for g in _group(np.array(yt_raw), gap=1) if 10 <= len(g) <= 30]
    yt_pos = sorted((g[0] + g[-1]) / 2.0 + 0.5 for g in yt_groups)
    y_ticks = D.regular_ticks(yt_pos, tol=4.0)
    if len(y_ticks) != 5:
        raise D.CalibrationError(f"Fig. 4B top phi_Lo: expected 5 y-ticks (0.50..0.30), got {y_ticks}")
    ya = D.Axis("linear", list(zip(y_ticks, [0.50, 0.45, 0.40, 0.35, 0.30])))

    return xa, ya, black, (c_tl0, c_tl1, r_top_x0)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, black, (c_tl0, c_tl1, r_top_x0) = _calibrate(image)
    H, W = image.shape[:2]

    # Restrict to top panel interior (exclude right of vertical dashed line at t ~ 13 min)
    r0_roi, r1_roi = 150, r_top_x0 - 4
    c0_roi, c1_roi = max(0, c_tl0 - 70), int(round(xa.to_pixel(12.2)))
    roi = black[r0_roi:r1_roi, c0_roi:c1_roi]

    # The black oval markers are ~100 px wide x ~55 px tall, whereas spines/ticks/lines are <= 19 px thick.
    # Separable 1D box erosion (29 x 55 px) isolates ONLY the symmetric cores of the 9 filled oval markers!
    er = ndimage.minimum_filter1d(roi, size=29, axis=0)
    er = ndimage.minimum_filter1d(er, size=55, axis=1)

    lab, n_lab = ndimage.label(er)
    pts = []
    for i in range(1, n_lab + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 100:
            continue
        r_mid = float(rr.mean()) + 0.5 + r0_roi
        c_mid = float(cc.mean()) + 0.5 + c0_roi
        pts.append((c_mid, r_mid))
    pts.sort(key=lambda p: p[0])

    if len(pts) != 9:
        raise RuntimeError(f"Fig. 4B phi_Lo: expected 9 markers, found {len(pts)}")

    err_px = 2.0 + max(ya.residual_px, 0.5)
    out = []
    for c_mid, r_mid in pts:
        x_val = float(xa.to_value(c_mid))
        y_val = float(ya.to_value(r_mid))
        a, b = ya.to_value(r_mid - err_px), ya.to_value(r_mid + err_px)
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [c_mid, r_mid],
            "px_lo": float(r_mid - err_px),
            "px_hi": float(r_mid + err_px),
            "note": f"oval marker; cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px",
        })
    return out
