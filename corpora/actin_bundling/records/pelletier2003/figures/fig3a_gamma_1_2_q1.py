"""Calibrated reader for Fig. 3(a) (p. 3): resolved first-order SAXS peak wavevector q_1 for gamma = 1:2."""

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

    xt = [x + int(c0) - 10 for x in _cluster(np.where(gray[int(r1) + 10 : int(r1) + 22, int(c0) - 10 : int(c1) + 15].mean(axis=0) < 180)[0])]
    xv = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0]
    ax_x = Axis("log10", list(zip(xt, xv)))
    return ax_x, (c0, r0, c1, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, (c0, r0, c1, r1) = _calibrate_fig3a(gray)

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
    mask = lbl == comps[0][1]  # gamma = 1:2 (top curve)

    # Trace curve centroid strictly right of the low-q valley (~0.165) across [0.172, 0.225]
    c_start = int(ax_x.to_pixel(0.172))
    c_end = int(ax_x.to_pixel(0.225))
    cols, rows = [], []
    for c in range(c_start, c_end + 1):
        rr = np.where(mask[:, c])[0]
        if len(rr):
            w = 255.0 - gray[rr, c]
            cols.append(float(c))
            rows.append(float(np.sum(rr * w) / np.sum(w)))
    cols_arr = np.array(cols)
    rows_arr = np.array(rows)
    r_min = float(rows_arr.min())
    crest = rows_arr <= r_min + 4.0  # within 1 native pixel of local maximum
    c_peak = float(cols_arr[crest].mean())
    r_peak = float(rows_arr[crest].mean())
    q_peak = float(ax_x.to_value(c_peak))

    px_quant = 2.0
    crest_half = 0.5 * float(cols_arr[crest].max() - cols_arr[crest].min())
    c_err = float(np.hypot(np.hypot(px_quant, ax_x.residual_px), max(2.0, 0.5 * crest_half)))
    q_lo = float(ax_x.to_value(c_peak - c_err))
    q_hi = float(ax_x.to_value(c_peak + c_err))

    return [{
        "x": float(at[0]) if at else 1.0,
        "y": q_peak,
        "lo": min(q_lo, q_hi),
        "hi": max(q_lo, q_hi),
        "px": [c_peak, r_peak],
        "note": f"raw local peak q_1 (gamma=1:2, unsubtracted); x_res={ax_x.residual_px:.2f}px, crest_w={2*crest_half:.0f}px",
    }]
