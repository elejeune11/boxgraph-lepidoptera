"""Calibrated reader for Fig. 2A (p. 2): 0 mM Ba2+ SAXS profile I(q)/I_0 vs q."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from boxgraph.digitize import Axis


def _longest_run_len(arr: np.ndarray) -> int:
    if not arr.any():
        return 0
    padded = np.concatenate(([False], arr, [False]))
    diff = np.diff(padded.astype(int))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    return int((ends - starts).max())


def _cluster_1d(idxs: np.ndarray, gap: int = 4) -> list[float]:
    if len(idxs) == 0:
        return []
    groups = np.split(idxs, np.where(np.diff(idxs) > gap)[0] + 1)
    return [float(g.mean()) for g in groups if len(g) > 0]


def _calibrate_fig2a(gray: np.ndarray):
    H, W = gray.shape
    dark = gray < 170
    # Solid frame lines span >70% of image width/height continuously
    h_rows = np.array([r for r in range(H) if _longest_run_len(dark[r]) > 0.70 * W])
    v_cols = np.array([c for c in range(W) if _longest_run_len(dark[:, c]) > 0.70 * H])
    hs = _cluster_1d(h_rows)
    vs = _cluster_1d(v_cols)
    r0, r1 = hs[0], hs[-1]
    c0, c1 = vs[0], vs[-1]

    # Bottom outward x-ticks (14 ticks: 0.05 to 0.18 in steps of 0.01 A^-1)
    bt = gray[int(r1) + 10 : int(r1) + 34, max(0, int(c0) - 12) : min(W, int(c1) + 12)] < 170
    xt = [x + max(0, int(c0) - 12) for x in _cluster_1d(np.where(bt.mean(axis=0) > 0.30)[0])]
    x_vals = [0.05 + 0.01 * i for i in range(len(xt))]
    ax_x = Axis("linear", list(zip(xt, x_vals)))

    # Left outward logarithmic y-ticks (0.6..0.9, 1..10 in relative decade units; lower major tick = 1)
    lt = gray[max(0, int(r0) - 10) : min(H, int(r1) + 10), max(0, int(c0) - 35) : max(1, int(c0) - 10)] < 175
    yt = [y + max(0, int(r0) - 10) for y in _cluster_1d(np.where(lt.mean(axis=1) > 0.25)[0])]
    # Identify the two long major decade ticks (len > 30 px left of c0)
    majors = []
    for y in yt:
        ry = int(round(y))
        strip = gray[ry - 2 : ry + 3, max(0, int(c0) - 60) : int(c0) - 6].mean(axis=0) < 175
        if int(strip.sum()) > 30:
            majors.append(y)
    # Match every detected tick to its nearest log10 minor/major value (10 at top major, 1 at bottom major)
    y_hi, y_lo = majors[0], majors[-1]
    px_per_dec = y_lo - y_hi
    cand_vals = [0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    matched = []
    for y in yt:
        est = 10.0 ** ((y_lo - y) / px_per_dec)
        best_v = min(cand_vals, key=lambda v: abs(np.log10(v) - np.log10(est)))
        if abs(y - (y_lo - px_per_dec * np.log10(best_v))) < 12.0:
            matched.append((y, best_v))
    ax_y = Axis("log10", matched)
    return ax_x, ax_y, (c0, r0, c1, r1)


def _extract_curve_mask(gray: np.ndarray, box: tuple[float, float, float, float], which_idx: int) -> np.ndarray:
    """which_idx: 0 = top (108 mM), 1 = middle (36 mM), 2 = bottom (0 mM)."""
    c0, r0, c1, r1 = box
    interior = np.zeros_like(gray, dtype=bool)
    rs, re = int(r0 + 35), int(r1 - 15)
    cs, ce = int(c0 + 10), int(c1 - 10)
    interior[rs:re, cs:ce] = gray[rs:re, cs:ce] < 165
    lbl, n = ndimage.label(interior)
    spanning = []
    for k in range(1, n + 1):
        rr, cc = np.where(lbl == k)
        if cc.max() - cc.min() > 0.70 * (c1 - c0):
            spanning.append((float(rr.mean()), k))
    spanning.sort()
    return lbl == spanning[which_idx][1]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, ax_y, box = _calibrate_fig2a(gray)
    mask = _extract_curve_mask(gray, box, which_idx=2)  # 0 mM (bottom curve)

    # Native 150-ppi raster upsampled 8x at 1200 dpi -> 0.5 native px = 4.0 px
    px_quant = 4.0
    x_err_px = np.hypot(px_quant, ax_x.residual_px)
    y_cal_px = np.hypot(px_quant, ax_y.residual_px)

    out = []
    for q in at:
        c = ax_x.to_pixel(q)
        ci = int(round(c))
        # Average over a +/- 4 px horizontal strip (~1 native pixel) for robust centroid
        rr_list, ww_list = [], []
        for dc in range(-4, 5):
            cc = ci + dc
            if 0 <= cc < gray.shape[1]:
                rr = np.where(mask[:, cc])[0]
                if len(rr):
                    rr_list.append(rr)
                    ww_list.append(255.0 - gray[rr, cc])
        if not rr_list:
            out.append({"x": float(q), "y": None, "note": "0 mM curve not found at q"})
            continue
        rr_all = np.concatenate(rr_list)
        ww_all = np.concatenate(ww_list)
        row = float(np.sum(rr_all * ww_all) / np.sum(ww_all))
        half_thick = 0.5 * float(rr_all.max() - rr_all.min())

        # Local slope estimate over +/- 12 px
        r_left = np.where(mask[:, max(0, ci - 12)])[0]
        r_right = np.where(mask[:, min(gray.shape[1] - 1, ci + 12)])[0]
        slope_px = 0.0
        if len(r_left) and len(r_right):
            slope_px = abs(float(r_right.mean() - r_left.mean())) / 24.0

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
