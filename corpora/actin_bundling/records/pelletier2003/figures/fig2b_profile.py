"""Calibrated reader for Fig. 2(b) (p. 3): LSCM fluorescence intensity line profile across crossing point b (gamma = 1:5)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from boxgraph.digitize import Axis


def _cluster(idxs: np.ndarray, gap: int = 4) -> list[float]:
    if len(idxs) == 0:
        return []
    groups = np.split(idxs, np.where(np.diff(idxs) > gap)[0] + 1)
    return [float(g.mean()) for g in groups if len(g) > 0]


def _calibrate_fig2b(img: np.ndarray):
    g = img.astype(float).mean(axis=2)
    dark = g < 150
    H, W = g.shape
    hs = _cluster(np.where(dark[:, 200 : W - 30].sum(axis=1) > 900)[0], gap=6)
    vs = _cluster(np.where(dark[50 : H - 150, :].sum(axis=0) > 900)[0], gap=6)
    r0, r1 = hs[0], hs[-1]
    c0, c1 = vs[0], vs[-1]

    # Y-axis labels: 250, 150, 50 (below "Int.")
    ylbl = (g[int(r0) - 20 : int(r1) + 20, max(0, int(c0) - 240) : int(c0) - 20] < 180).sum(axis=1)
    yrows = np.where(ylbl > 3)[0] + (int(r0) - 20)
    ygrs = [float(gr.mean()) for gr in np.split(yrows, np.where(np.diff(yrows) > 35)[0] + 1) if len(gr) > 15]
    ax_y_lbl = Axis("linear", list(zip(ygrs[-3:], [250.0, 150.0, 50.0])))

    lt = g[int(r0) - 4 : int(r1) + 4, int(c0) - 14 : int(c0) - 7].mean(axis=1)
    yt_cands = [y + int(r0) - 4 for y in _cluster(np.where(lt < 232)[0], gap=4)]
    ref_y = []
    for v in range(0, 301, 50):
        pred = ax_y_lbl.to_pixel(v)
        near = [t for t in yt_cands if abs(t - pred) < 16]
        if near:
            ref_y.append((min(near, key=lambda t: abs(t - pred)), float(v)))
    ax_y = Axis("linear", ref_y) if len(ref_y) >= 3 else ax_y_lbl

    # X-axis labels: 0.0, 0.4, 0.8, 1.2
    xlbl = (g[int(r1) + 25 : int(r1) + 110, :] < 180).sum(axis=0)
    xcols = np.where(xlbl > 3)[0]
    xgrs = [float(gr.mean()) for gr in np.split(xcols, np.where(np.diff(xcols) > 65)[0] + 1) if len(gr) > 30]
    ax_x_lbl = Axis("linear", list(zip(xgrs[:4], [0.0, 0.4, 0.8, 1.2])))

    bt = g[int(r1) + 10 : int(r1) + 15, int(c0) - 4 : int(c1) + 5].mean(axis=0)
    xt_cands = [x + int(c0) - 4 for x in _cluster(np.where(bt < 232)[0], gap=4)]
    ref_x = [(c0, 0.0)]
    for v in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]:
        pred = ax_x_lbl.to_pixel(v)
        near = [t for t in xt_cands if abs(t - pred) < 16]
        if near:
            ref_x.append((min(near, key=lambda t: abs(t - pred)), v))
    ax_x = Axis("linear", ref_x) if len(ref_x) >= 3 else ax_x_lbl
    return ax_x, ax_y, (c0, r0, c1, r1)


def _extract_profile_mask(img: np.ndarray, box: tuple[float, float, float, float]) -> np.ndarray:
    c0, r0, c1, r1 = box
    g = img.astype(float).mean(axis=2)
    mask = np.zeros_like(g, dtype=bool)
    rs, re = int(r0 + 12), int(r1 - 12)
    cs, ce = int(c0 + 12), int(c1 - 12)
    mask[rs:re, cs:ce] = g[rs:re, cs:ce] < 165

    # Remove the top-right inset micrograph precisely (tight 4px pad to avoid curve near x=0.80)
    green = (img[:, :, 1].astype(float) - img[:, :, 0].astype(float)) > 25
    lbl, _ = ndimage.label(mask)
    inset_labels = set(np.unique(lbl[green])) - {0}
    for k in inset_labels:
        rr, cc = np.where(lbl == k)
        mask[max(rs, rr.min() - 4) : min(re, rr.max() + 5), max(cs, cc.min() - 4) : min(ce, cc.max() + 5)] = False
    # Keep only the continuous line-profile component (rejects any isolated stray specks)
    lbl2, n2 = ndimage.label(mask)
    if n2 > 0:
        sizes = ndimage.sum(mask, lbl2, range(1, n2 + 1))
        mask = lbl2 == (int(np.argmax(sizes)) + 1)
    return mask


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, ax_y, box = _calibrate_fig2b(image)
    mask = _extract_profile_mask(image, box)

    px_quant = 4.0
    x_err_px = float(np.hypot(px_quant, ax_x.residual_px))
    y_cal_px = float(np.hypot(px_quant, ax_y.residual_px))

    out = []
    for xv in at:
        c = ax_x.to_pixel(xv)
        ci = int(round(c))
        rr_list, ww_list = [], []
        for dc in range(-4, 5):
            cc = ci + dc
            if 0 <= cc < gray.shape[1]:
                rr = np.where(mask[:, cc])[0]
                if len(rr):
                    rr_list.append(rr)
                    ww_list.append(255.0 - gray[rr, cc])
        if not rr_list:
            out.append({"x": float(xv), "y": None, "note": "profile not found at x"})
            continue
        rr_all = np.concatenate(rr_list)
        ww_all = np.concatenate(ww_list)
        row = float(np.sum(rr_all * ww_all) / np.sum(ww_all))
        half_thick = 0.5 * float(rr_all.max() - rr_all.min())

        r_left = np.where(mask[:, max(0, ci - 16)])[0]
        r_right = np.where(mask[:, min(gray.shape[1] - 1, ci + 16)])[0]
        slope_px = abs(float(r_right.mean() - r_left.mean())) / 32.0 if (len(r_left) and len(r_right)) else 0.0

        row_err = float(np.hypot(np.hypot(y_cal_px, 0.5 * half_thick), slope_px * x_err_px))
        y_val = ax_y.to_value(row)
        y_hi = ax_y.to_value(row - row_err)
        y_lo = ax_y.to_value(row + row_err)
        out.append({
            "x": float(xv),
            "y": float(y_val),
            "lo": float(min(y_lo, y_hi)),
            "hi": float(max(y_lo, y_hi)),
            "px": [float(c), float(row)],
            "px_lo": float(row - row_err),
            "px_hi": float(row + row_err),
            "note": f"x_res={ax_x.residual_px:.2f}px, y_res={ax_y.residual_px:.2f}px",
        })
    return out
