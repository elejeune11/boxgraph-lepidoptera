"""Calibrated reader for Fig. 3(b) (p. 3): plotted 1/q rod-scattering background fit curve (gamma = 1:5)."""

from __future__ import annotations

import numpy as np

from boxgraph.digitize import Axis


def _cluster(idxs: np.ndarray, gap: int = 4) -> list[float]:
    if len(idxs) == 0:
        return []
    groups = np.split(idxs, np.where(np.diff(idxs) > gap)[0] + 1)
    return [float(g.mean()) for g in groups if len(g) > 0]


def _calibrate_fig3b(gray: np.ndarray):
    dark = gray < 140
    r1 = _cluster(np.where(dark[100:1700, 200:950].sum(axis=1) > 500)[0])[-1] + 100.0
    c0 = _cluster(np.where(dark[200:1550, 100:900].sum(axis=0) > 800)[0])[0] + 100.0

    xt_out = [x + int(c0) + 30 for x in _cluster(np.where(gray[int(r1) + 7 : int(r1) + 10, int(c0) + 30 : 950].mean(axis=0) < 200)[0])]
    xt = [c0] + xt_out[:3]
    ax_x = Axis("linear", list(zip(xt, [0.1, 0.2, 0.3, 0.4])))

    yt_maj = [y + 200.0 for y in _cluster(np.where(gray[200 : int(r1) + 10, int(c0) - 22 : int(c0) - 15].mean(axis=1) < 200)[0])]
    ax_y = Axis("linear", list(zip(yt_maj, [3000.0, 2000.0, 1000.0])))
    return ax_x, ax_y, (c0, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    gray = image.astype(float).mean(axis=2)
    ax_x, ax_y, (c0, r1) = _calibrate_fig3b(gray)

    px_quant = 2.0
    x_err_px = float(np.hypot(px_quant, ax_x.residual_px))
    y_cal_px = float(np.hypot(px_quant, ax_y.residual_px))

    out = []
    for q in at:
        c = ax_x.to_pixel(q)
        ci = int(round(c))
        r_top = 600 if ci > 480 else 400
        r_bot = int(r1) - 10
        col = gray[r_top:r_bot, max(0, ci - 2) : min(gray.shape[1], ci + 3)].mean(axis=1)
        rr = np.where(col < 165)[0] + r_top
        if not len(rr):
            out.append({"x": float(q), "y": None, "note": "1/q background curve not found at q"})
            continue
        segs = np.split(rr, np.where(np.diff(rr) > 6)[0] + 1)
        if len(segs) < 2:
            out.append({"x": float(q), "y": None, "note": "background curve merged with data curve at q"})
            continue
        # Lower thin smooth 1/q fit is the bottommost segment
        seg = segs[-1]
        w = 255.0 - col[seg - r_top]
        row = float(np.sum(seg * w) / np.sum(w))
        half_thick = 0.5 * float(seg.max() - seg.min())

        slope_px = 0.0
        r_lr = []
        for cc_s in (ci - 12, ci + 12):
            rt_s = 600 if cc_s > 480 else 400
            if int(c0) + 6 < cc_s < gray.shape[1] - 4:
                rr_s = np.where(gray[rt_s:r_bot, cc_s] < 165)[0] + rt_s
                if len(rr_s):
                    ss = np.split(rr_s, np.where(np.diff(rr_s) > 6)[0] + 1)
                    if len(ss) >= 2:
                        r_lr.append(float(ss[-1].mean()))
        if len(r_lr) == 2:
            slope_px = abs(r_lr[1] - r_lr[0]) / 24.0

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
