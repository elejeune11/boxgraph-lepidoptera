"""Calibrated reader for Fig. 3(c) (p. 3): plotted second-order Gaussian component fit G_2(q) (gamma = 1:5)."""

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
        if 0.150 <= q <= 0.218:
            # Left of G1/G2 crossing: G_2 is the bottommost thin curve cleanly above x-axis (rows 1110..r_bot)
            r_lo_s, r_hi_s = 1110, r_bot
            col = gray[r_lo_s:r_hi_s, max(0, ci - 1) : min(gray.shape[1], ci + 2)].mean(axis=1)
            rr = np.where(col < 175)[0] + r_lo_s
            segs = np.split(rr, np.where(np.diff(rr) > 6)[0] + 1) if len(rr) else []
            if not segs:
                out.append({"x": float(q), "y": None, "note": "G2 not found on left wing"})
                continue
            seg = segs[-1]
        elif 0.236 <= q <= 0.265:
            # Between G1/G2 crossing and G2 peak: G_2 is the middle curve (above G1 tail, below data/sum)
            r_lo_s, r_hi_s = 700, r_bot
            col = gray[r_lo_s:r_hi_s, max(0, ci - 1) : min(gray.shape[1], ci + 2)].mean(axis=1)
            rr = np.where(col < 170)[0] + r_lo_s
            segs = np.split(rr, np.where(np.diff(rr) > 6)[0] + 1) if len(rr) else []
            if len(segs) < 3:
                out.append({"x": float(q), "y": None, "note": "G2 not cleanly isolated between data and G1"})
                continue
            seg = segs[1]
        elif 0.280 <= q <= 0.455:
            # Peak and right flank (G1 < 1 here, so G_2 coincides with the plotted sum/data curve)
            r_lo_s, r_hi_s = 600, r_bot
            col = gray[r_lo_s:r_hi_s, max(0, ci - 2) : min(gray.shape[1], ci + 3)].mean(axis=1)
            rr = np.where(col < 160)[0] + r_lo_s
            segs = np.split(rr, np.where(np.diff(rr) > 12)[0] + 1) if len(rr) else []
            if not segs:
                out.append({"x": float(q), "y": None, "note": "G2 not found on right flank"})
                continue
            seg = segs[0]
        else:
            out.append({"x": float(q), "y": None, "note": "G2 crosses G1 or merges near shoulder at this q"})
            continue

        w = 255.0 - col[seg - r_lo_s]
        row = float(np.sum(seg * w) / np.sum(w))
        half_thick = 0.5 * float(seg.max() - seg.min())

        # Local slope from +/-8 px
        slope_px = 0.0
        r_lr = []
        for cc_s in (ci - 8, ci + 8):
            if int(c0) + 6 < cc_s < gray.shape[1] - 4:
                rr_s = np.where(gray[r_lo_s:r_hi_s, cc_s] < 170)[0] + r_lo_s
                if len(rr_s):
                    ss = np.split(rr_s, np.where(np.diff(rr_s) > 6)[0] + 1)
                    if 0.145 <= q <= 0.218:
                        r_lr.append(float(ss[-1].mean()))
                    elif 0.236 <= q <= 0.265 and len(ss) >= 3:
                        r_lr.append(float(ss[1].mean()))
                    elif q >= 0.280:
                        r_lr.append(float(ss[0].mean()))
        if len(r_lr) == 2:
            slope_px = abs(r_lr[1] - r_lr[0]) / 16.0

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
