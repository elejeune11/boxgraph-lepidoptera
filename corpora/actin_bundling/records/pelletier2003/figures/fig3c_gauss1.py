"""Calibrated reader for Fig. 3(c) (p. 3): plotted first-order Gaussian component fit G_1(q) (gamma = 1:5)."""

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

    xt = [x + int(c0) - 15 for x in _cluster(np.where(gray[int(r1) + 8 : int(r1) + 18, int(c0) - 15 : 950].mean(axis=0) < 200)[0])]
    ax_x = Axis("linear", list(zip(xt[:4], [0.1, 0.2, 0.3, 0.4])))

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
        r_bot = int(r1) - 18
        # Isolate G_1 cleanly by region relative to data/sum (top) and G_2 (crossing near q=0.226)
        if 0.178 <= q <= 0.218:
            # G_1 is the middle curve below data dots (row > 555) and above G_2 (row < 1100)
            col = gray[555:1100, max(0, ci - 1) : min(gray.shape[1], ci + 2)].mean(axis=1)
            rr = np.where(col < 170)[0] + 555
            segs = np.split(rr, np.where(np.diff(rr) > 6)[0] + 1) if len(rr) else []
            if not segs:
                out.append({"x": float(q), "y": None, "note": "G1 not found in peak region"})
                continue
            seg = segs[-1]  # lowest segment in [555..1100] avoids any upper data-dot tail near 560
            r_ref = 555
        elif 0.232 <= q <= 0.258:
            # Right of G1/G2 crossing: G_1 is the only thin curve in [1130..r_bot] (above x-axis)
            col = gray[1130:r_bot, max(0, ci - 1) : min(gray.shape[1], ci + 2)].mean(axis=1)
            rr = np.where(col < 175)[0] + 1130
            segs = np.split(rr, np.where(np.diff(rr) > 6)[0] + 1) if len(rr) else []
            if not segs:
                out.append({"x": float(q), "y": None, "note": "G1 not found on right flank"})
                continue
            seg = segs[-1]
            r_ref = 1130
        else:
            out.append({"x": float(q), "y": None, "note": "G1 overlaps data or G2 crossing at this q"})
            continue

        w = 255.0 - col[seg - r_ref]
        row = float(np.sum(seg * w) / np.sum(w))
        half_thick = 0.5 * float(seg.max() - seg.min())

        # Estimate local slope numerically from +/-6 px
        slope_px = 0.0
        r_lr = []
        for cc_s in (ci - 6, ci + 6):
            r_lo_s, r_hi_s = (555, 1100) if q <= 0.218 else (1130, r_bot)
            rr_s = np.where(gray[r_lo_s:r_hi_s, cc_s] < 175)[0] + r_lo_s
            if len(rr_s):
                r_lr.append(float(np.split(rr_s, np.where(np.diff(rr_s) > 6)[0] + 1)[-1].mean()))
        if len(r_lr) == 2:
            slope_px = abs(r_lr[1] - r_lr[0]) / 12.0

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
