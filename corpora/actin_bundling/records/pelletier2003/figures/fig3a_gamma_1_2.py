"""Calibrated reader for Fig. 3(a) (p. 3): log-log SAXS I(q) vs q for gamma = 1:2 (top curve)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from boxgraph.digitize import Axis


def _cluster(idxs: np.ndarray, gap: int = 4) -> list[float]:
    if len(idxs) == 0:
        return []
    groups = np.split(idxs, np.where(np.diff(idxs) > gap)[0] + 1)
    return [float(g.mean()) for g in groups if len(g) > 0]


def _calibrate_fig3a(gray: np.ndarray):
    H, W = gray.shape
    dark = gray < 140
    hs = _cluster(np.where(dark[:, 100 : W - 50].sum(axis=1) > 1200)[0])
    vs = _cluster(np.where(dark[100 : H - 200, :].sum(axis=0) > 2500)[0])
    r0, r1 = hs[0], hs[-1]
    c0, c1 = vs[0], vs[-1]

    # Bottom outward ticks: 0.05..0.09, 0.1..0.9, 1..5 (19 logarithmic ticks)
    xt = [x + int(c0) - 10 for x in _cluster(np.where(gray[int(r1) + 10 : int(r1) + 22, int(c0) - 10 : int(c1) + 15].mean(axis=0) < 180)[0])]
    xv = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0]
    ax_x = Axis("log10", list(zip(xt, xv)))

    # Left outward major ticks: 10^4, 10^3, 10^2, 10^1, 10^0, 10^-1 (top border is 9e4)
    yt_maj = [y + int(r0) - 10 for y in _cluster(np.where(gray[int(r0) - 10 : int(r1) + 10, int(c0) - 35 : int(c0) - 25].mean(axis=1) < 180)[0])]
    ax_y0 = Axis("log10", list(zip(yt_maj, [1e4, 1e3, 1e2, 1e1, 1e0, 1e-1])))
    yt_all = [y + int(r0) - 10 for y in _cluster(np.where(gray[int(r0) - 10 : int(r1) + 10, int(c0) - 22 : int(c0) - 10].mean(axis=1) < 180)[0])]
    cand_y = [m * (10.0**e) for e in range(-1, 5) for m in range(1, 10)] + [1e5]
    matched_y = []
    for y in yt_all:
        est = ax_y0.to_value(y)
        best = min(cand_y, key=lambda v: abs(np.log10(v) - np.log10(est)))
        if abs(y - ax_y0.to_pixel(best)) < 10.0:
            matched_y.append((y, best))
    ax_y = Axis("log10", matched_y)
    return ax_x, ax_y, (c0, r0, c1, r1)


def _extract_curve_mask(gray: np.ndarray, box: tuple[float, float, float, float], which_idx: int) -> np.ndarray:
    c0, r0, c1, r1 = box
    interior = np.zeros_like(gray, dtype=bool)
    interior[int(r0) + 25 : int(r1) - 25, int(c0) + 20 : int(c1) - 35] = (
        gray[int(r0) + 25 : int(r1) - 25, int(c0) + 20 : int(c1) - 35] < 160
    )
    lbl, n = ndimage.label(interior)
    comps = []
    for k in range(1, n + 1):
        rr, cc = np.where(lbl == k)
        if len(rr) > 5000 and (cc.max() - cc.min()) > 600:
            comps.append((float(rr.mean()), k))
    comps.sort()
    return lbl == comps[which_idx][1]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, ax_y, box = _calibrate_fig3a(gray)
    mask = _extract_curve_mask(gray, box, which_idx=0)  # gamma = 1:2 (top curve)

    px_quant = 2.0  # 0.5 native pixel at 300 ppi -> 1200 dpi (4x)
    x_err_px = float(np.hypot(px_quant, ax_x.residual_px))
    y_cal_px = float(np.hypot(px_quant, ax_y.residual_px))

    out = []
    for q in at:
        c = ax_x.to_pixel(q)
        ci = int(round(c))
        rr_list, ww_list = [], []
        for dc in range(-3, 4):
            cc = ci + dc
            if 0 <= cc < gray.shape[1]:
                rr = np.where(mask[:, cc])[0]
                if len(rr):
                    rr_list.append(rr)
                    ww_list.append(255.0 - gray[rr, cc])
        if not rr_list:
            out.append({"x": float(q), "y": None, "note": "gamma=1:2 curve not found at q"})
            continue
        rr_all = np.concatenate(rr_list)
        ww_all = np.concatenate(ww_list)
        row = float(np.sum(rr_all * ww_all) / np.sum(ww_all))
        half_thick = 0.5 * float(rr_all.max() - rr_all.min())

        r_left = np.where(mask[:, max(0, ci - 10)])[0]
        r_right = np.where(mask[:, min(gray.shape[1] - 1, ci + 10)])[0]
        slope_px = abs(float(r_right.mean() - r_left.mean())) / 20.0 if (len(r_left) and len(r_right)) else 0.0

        row_err = float(np.hypot(np.hypot(y_cal_px, 0.5 * half_thick), slope_px * x_err_px))
        y_val = ax_y.to_value(row)
        y_hi = ax_y.to_value(row - row_err)
        y_lo = ax_y.to_value(row + row_err)
        out.append({
            "x": float(q),
            "y": float(y_val),
            "lo": float(min(y_lo, y_hi)),
            "hi": float(max(y_lo, y_hi)),
            "px": [float(c), float(row)],
            "px_lo": float(row - row_err),
            "px_hi": float(row + row_err),
            "note": f"x_res={ax_x.residual_px:.2f}px, y_res={ax_y.residual_px:.2f}px",
        })
    return out
