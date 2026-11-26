"""Calibrated reader for Fig. 3(c) (p. 3): background-subtracted SAXS data I_sub(q) for gamma = 1:5."""

from __future__ import annotations

import numpy as np

from boxgraph.digitize import Axis


def _cluster(idxs: np.ndarray, gap: int = 4) -> list[float]:
    if len(idxs) == 0:
        return []
    groups = np.split(idxs, np.where(np.diff(idxs) > gap)[0] + 1)
    return [float(g.mean()) for g in groups if len(g) > 0]


def _calibrate_fig3c(gray: np.ndarray):
    dark = gray < 140
    r1 = _cluster(np.where(dark[100:1700, 200:950].sum(axis=1) > 500)[0])[-1] + 100.0
    c0 = _cluster(np.where(dark[200:1550, 100:900].sum(axis=0) > 800)[0])[0] + 100.0

    # Bottom outward ticks: 0.1, 0.2, 0.3, 0.4
    xt = [x + int(c0) - 15 for x in _cluster(np.where(gray[int(r1) + 8 : int(r1) + 18, int(c0) - 15 : 950].mean(axis=0) < 200)[0])]
    ax_x = Axis("linear", list(zip(xt[:4], [0.1, 0.2, 0.3, 0.4])))

    # Left outward ticks: 8, 6, 4, 2, 0 (x 10^2 -> 800, 600, 400, 200, 0)
    yt = [y + 200.0 for y in _cluster(np.where(gray[200 : int(r1) + 10, int(c0) - 16 : int(c0) - 6].mean(axis=1) < 200)[0])]
    ax_y = Axis("linear", list(zip(yt[:5], [800.0, 600.0, 400.0, 200.0, 0.0])))
    return ax_x, ax_y, (c0, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, ax_y, (c0, r1) = _calibrate_fig3c(gray)

    px_quant = 2.0
    x_err_px = float(np.hypot(px_quant, ax_x.residual_px))
    y_cal_px = float(np.hypot(px_quant, ax_y.residual_px))

    out = []
    for q in at:
        c = ax_x.to_pixel(q)
        ci = int(round(c))
        # Mask out panel label "(c)" at rows < 560 for c > 620
        r_top = 560 if ci > 620 else 280
        r_bot = int(r1) - 6
        # Average across narrow 3-px column band so discrete markers connect cleanly
        col = gray[r_top:r_bot, max(0, ci - 2) : min(gray.shape[1], ci + 3)].mean(axis=1)
        rr = np.where(col < 160)[0] + r_top
        if not len(rr):
            out.append({"x": float(q), "y": None, "note": "background-subtracted curve not found at q"})
            continue
        # Merge closely adjacent data-dot clusters within 30 px on the steep rising/falling flanks,
        # while keeping distinct lower Gaussian component curves separated
        segs = np.split(rr, np.where(np.diff(rr) > 22)[0] + 1)
        seg = segs[0]
        # On the very steep rising flank (0.145 <= q <= 0.175), markers are staggered vertically by ~25-40 px;
        # restrict to within 40 px of the topmost dark ink to avoid including G1 below
        if len(segs) > 1 and (segs[1][0] - seg[-1]) < 45 and 0.140 <= q <= 0.175:
            seg = np.concatenate([seg, segs[1]])
        w = 255.0 - col[seg - r_top]
        row = float(np.sum(seg * w) / np.sum(w))
        half_thick = 0.5 * float(seg.max() - seg.min())

        slope_px = 0.0
        r_lr = []
        for cc_s in (ci - 10, ci + 10):
            rt_s = 560 if cc_s > 620 else 280
            if int(c0) + 6 < cc_s < gray.shape[1] - 4:
                c_s = gray[rt_s:r_bot, max(0, cc_s - 2) : min(gray.shape[1], cc_s + 3)].mean(axis=1)
                rr_s = np.where(c_s < 160)[0] + rt_s
                if len(rr_s):
                    s0 = np.split(rr_s, np.where(np.diff(rr_s) > 22)[0] + 1)[0]
                    r_lr.append(float(s0.mean()))
        if len(r_lr) == 2:
            slope_px = abs(r_lr[1] - r_lr[0]) / 20.0

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
