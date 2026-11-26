"""Calibrated point reader for si_fig1a_D_vs_t (SI Fig. 1, p. 20)."""

import math
import numpy as np
from scipy import ndimage

Y_TICK_VALS = [0.1, 0.08, 0.06, 0.04, 0.02, 0.0, -0.02]


def _groups(arr):
    gs = [[int(arr[0])]]
    for v in arr[1:]:
        iv = int(v)
        if iv == gs[-1][-1] + 1:
            gs[-1].append(iv)
        else:
            gs.append([iv])
    return [(g[0], g[-1], float(np.mean(g)) + 0.5) for g in gs]


def _interior_ticks(strip_2d, offset):
    col_dark = strip_2d[:8, :].mean(axis=0) > 0.6
    idx = np.where(col_dark)[0]
    gs = [[int(idx[0])]]
    for v in idx[1:]:
        iv = int(v)
        if iv <= gs[-1][-1] + 2:
            gs[-1].append(iv)
        else:
            gs.append([iv])
    return [offset + float(np.mean(g)) + 0.5 for g in gs if len(g) <= 15]


def read(image, at):
    blue = (image[:, :, 2] > 150) & (image[:, :, 0] < 100) & (image[:, :, 1] < 100)
    hg = _groups(np.where(blue.mean(axis=1) > 0.5)[0])
    vg = _groups(np.where(blue.mean(axis=0) > 0.5)[0])
    r0_a, r0_b, r0 = hg[0]
    r1_a, r1_b, r1 = hg[-1]
    c0_a, c0_b, c0 = vg[0]
    c1_a, c1_b, c1 = vg[-1]

    x_ticks = [c0] + _interior_ticks(blue[r1_a - 35:r1_a, c0_b + 1:c1_a][::-1, :], c0_b + 1) + [c1]
    y_ticks = [r0] + _interior_ticks(blue[r0_b + 1:r1_a, c1_a - 35:c1_a][:, ::-1].T, r0_b + 1) + [r1]
    x_vals = list(range(len(x_ticks)))
    bx, ax = np.polyfit(x_ticks, x_vals, 1)
    by, ay = np.polyfit(y_ticks, Y_TICK_VALS, 1)
    x_rms = float(np.sqrt(np.mean((np.array(x_ticks) - (np.array(x_vals) - ax) / bx) ** 2)))
    y_rms = float(np.sqrt(np.mean((np.array(y_ticks) - (np.array(Y_TICK_VALS) - ay) / by) ** 2)))

    sub = image[r0_b + 3:r1_a - 3, c0_b + 3:c1_a - 3]
    red = (sub[:, :, 0] > 180) & (sub[:, :, 1] < 80) & (sub[:, :, 2] < 80)
    lbl, n = ndimage.label(red)
    pts = []
    for i in range(1, n + 1):
        rr, cc = np.where(lbl == i)
        if len(rr) < 20:
            continue
        cy = (r0_b + 3) + float(rr.mean()) + 0.5
        cx = (c0_b + 3) + float(cc.mean()) + 0.5
        pts.append((cx, cy, len(rr)))
    pts.sort(key=lambda p: p[0])

    half_px = math.sqrt(2.0 ** 2 + (2.0 * y_rms) ** 2)
    out = []
    for cx, cy, area in pts:
        xv = ax + bx * cx
        yv = ay + by * cy
        y_hi = ay + by * (cy - half_px)
        y_lo = ay + by * (cy + half_px)
        out.append({
            "x": float(xv),
            "y": float(yv),
            "lo": float(y_lo),
            "hi": float(y_hi),
            "px": [float(cx), float(cy)],
            "px_lo": float(cy - half_px),
            "px_hi": float(cy + half_px),
            "note": f"red point (area={area}); x_rms={x_rms:.2f}px, y_rms={y_rms:.2f}px",
        })
    return out
