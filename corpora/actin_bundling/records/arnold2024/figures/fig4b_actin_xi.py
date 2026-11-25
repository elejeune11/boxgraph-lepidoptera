"""Calibrated reader for Fig. 4B (top panel, right y-axis, red; PDF p. 6 / printed p. 26575):
Actin mean wavelength xi [um] vs. time t [min] during cholesterol insertion (t = 0 to 18 min).

Note on right y-axis ticks: the red ticks on the right spine show ~10-13 px (~0.002 um)
typesetting jitter from uniform spacing; fitting all 7 detected ticks (3.00 to 2.70 um in 0.05 um
steps) automatically propagates this residual into the 95% reading interval.
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
    red = (r > 160) & (g < 100) & (b < 100)

    # Horizontal black x-axes (top panel and bottom panel)
    hg = _group(np.where(black.sum(axis=1) > 0.5 * W)[0], gap=6)
    r_top_x0, r_top_x1 = hg[0][0], hg[0][-1]
    r_bot_x0, r_bot_x1 = hg[1][0], hg[1][-1]
    r_top_xmid = (r_top_x0 + r_top_x1) / 2.0 + 0.5

    # Left black spine & right red spine of top panel
    tl_g = _group(np.where(black[:r_top_x0, 600:1000].sum(axis=0) > 0.5 * r_top_x0)[0] + 600, gap=4)[0]
    tr_g = _group(np.where(red[:r_top_x0, 4000:5000].sum(axis=0) > 0.5 * r_top_x0)[0] + 4000, gap=4)[0]
    c_tl1 = tl_g[-1]
    c_zero = (tl_g[0] + tl_g[-1]) / 2.0 + 0.5
    c_tr0 = tr_g[0]

    # X-axis ticks from bottom panel (shared exact columns: 0, 5, 10, 15, 20 min)
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
        raise D.CalibrationError(f"Fig. 4B top actin_xi: expected 5 x-ticks, got {x_ticks}")
    xa = D.Axis("linear", list(zip(x_ticks, [0.0, 5.0, 10.0, 15.0, 20.0])))

    # Right red y-axis ticks (pointing leftward from right red spine: 3.00, 2.95, 2.90, 2.85, 2.80, 2.75)
    yt_raw = []
    for r_idx in range(100, r_top_x0 - 10):
        n = 0
        while c_tr0 - 1 - n >= 0 and red[r_idx, c_tr0 - 1 - n]:
            n += 1
        if 25 <= n <= 100:
            yt_raw.append(r_idx)
    yt_groups = [g for g in _group(np.array(yt_raw), gap=1) if 10 <= len(g) <= 30]
    yt_pos = sorted((g[0] + g[-1]) / 2.0 + 0.5 for g in yt_groups)
    if len(yt_pos) != 6:
        raise D.CalibrationError(f"Fig. 4B top actin_xi: expected 6 red y-ticks above bottom spine, got {yt_pos}")
    # Include bottom spine (2.70 um) as the 7th mark
    ya = D.Axis("linear", list(zip(yt_pos + [r_top_xmid], [3.00, 2.95, 2.90, 2.85, 2.80, 2.75, 2.70])))

    return xa, ya, red, (tl_g[0], c_tr0, r_top_x0)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, red, (c_tl0, c_tr0, r_top_x0) = _calibrate(image)
    H, W = image.shape[:2]

    # Red oval markers have lighter/hollow centers (~88-120 px wide x ~50-54 px tall);
    # fill marker interiors before separable box erosion (29 x 53 px).
    r0_roi, r1_roi = 150, r_top_x0 - 4
    c0_roi, c1_roi = max(0, c_tl0 - 70), c_tr0 - 80
    r_ch, g_ch, b_ch = image[r0_roi:r1_roi, c0_roi:c1_roi, 0], image[r0_roi:r1_roi, c0_roi:c1_roi, 1], image[r0_roi:r1_roi, c0_roi:c1_roi, 2]
    roi = (r_ch > 160) & (g_ch < 120) & (b_ch < 120)
    filled = ndimage.binary_fill_holes(roi)
    er = ndimage.minimum_filter1d(filled, size=29, axis=0)
    er = ndimage.minimum_filter1d(er, size=53, axis=1)

    lab, n_lab = ndimage.label(er)
    pts = []
    for i in range(1, n_lab + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 50:
            continue
        r_mid = float(rr.mean()) + 0.5 + r0_roi
        c_mid = float(cc.mean()) + 0.5 + c0_roi
        pts.append((c_mid, r_mid))
    pts.sort(key=lambda p: p[0])

    if len(pts) != 12:
        raise RuntimeError(f"Fig. 4B actin_xi: expected 12 red markers, found {len(pts)}")

    # Propagate the right red axis tick-placement residual (~12.8 px) into the 95% interval
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
            "note": f"red oval marker; cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px",
        })
    return out
