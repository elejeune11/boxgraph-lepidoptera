"""Calibrated orange-filled-square marker reader for fig3b_data_fascin_002 (Fig. 3B, p. 19)."""

import math
import numpy as np
from scipy import ndimage, signal


def _groups(arr):
    if len(arr) == 0:
        return []
    gs = [[int(arr[0])]]
    for v in arr[1:]:
        iv = int(v)
        if iv == gs[-1][-1] + 1:
            gs[-1].append(iv)
        else:
            gs.append([iv])
    return [(g[0], g[-1], float(np.mean(g)) + 0.5) for g in gs]


def _tick_positions(strip_2d, offset):
    col_dark = strip_2d[:8, :].mean(axis=0) > 0.6
    idx = np.where(col_dark)[0]
    if len(idx) == 0:
        return []
    gs = [[int(idx[0])]]
    for v in idx[1:]:
        iv = int(v)
        if iv <= gs[-1][-1] + 2:
            gs[-1].append(iv)
        else:
            gs.append([iv])
    out = []
    for g in gs:
        if len(g) > 15:
            continue
        sub = strip_2d[:, g[0]:g[-1] + 1].any(axis=1)
        depth = 0
        for v in sub:
            if v:
                depth += 1
            else:
                break
        if 10 <= depth <= 28:
            out.append(offset + float(np.mean(g)) + 0.5)
    return out


def _calibrate(image):
    dark = image.mean(axis=2) < 200
    hg = _groups(np.where(dark.mean(axis=1) > 0.5)[0])
    vg = _groups(np.where(dark.mean(axis=0) > 0.5)[0])
    r0_a, r0_b, r0 = hg[0]
    r1_a, r1_b, r1 = hg[-1]
    c0_a, c0_b, c0 = vg[0]
    c1_a, c1_b, c1 = vg[-1]

    bx_raw = _tick_positions(dark[r1_a - 35:r1_a, c0_b + 1:c1_a][::-1, :], c0_b + 1)
    b0 = math.log10(40.0) / (c1 - c0)
    a0 = -b0 * c0
    xm_p, xm_v = [c0, c1], [1.0, 40.0]
    for v in (2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30):
        p_exp = (math.log10(v) - a0) / b0
        near = [p for p in bx_raw if abs(p - p_exp) < 15]
        if near:
            xm_p.append(min(near, key=lambda p: abs(p - p_exp)))
            xm_v.append(float(v))
    bx, ax = np.polyfit(xm_p, np.log10(xm_v), 1)
    x_rms = float(np.sqrt(np.mean((np.array(xm_p) - (np.log10(xm_v) - ax) / bx) ** 2)))

    ly_raw = _tick_positions(dark[r0_b + 1:r1_a, c0_b + 1:c0_b + 35].T, r0_b + 1)
    ry_raw = _tick_positions(dark[r0_b + 1:r1_a, c1_a - 35:c1_a][:, ::-1].T, r0_b + 1)
    by0 = (math.log10(0.04) - math.log10(80.0)) / (r1 - r0)
    ay0 = math.log10(80.0) - by0 * r0
    y_cands = (
        [0.05, 0.06, 0.07, 0.08, 0.09]
        + [k * 0.1 for k in range(1, 10)]
        + list(range(1, 10))
        + [k * 10 for k in range(1, 8)]
    )
    ym_p, ym_v = [], []
    for raw in (ly_raw, ry_raw):
        for v in y_cands:
            p_exp = (math.log10(v) - ay0) / by0
            near = [p for p in raw if abs(p - p_exp) < 18]
            if near:
                ym_p.append(min(near, key=lambda p: abs(p - p_exp)))
                ym_v.append(float(v))
    by, ay = np.polyfit(ym_p, np.log10(ym_v), 1)
    y_rms = float(np.sqrt(np.mean((np.array(ym_p) - (np.log10(ym_v) - ay) / by) ** 2)))
    return (ax, bx, x_rms, len(xm_p)), (ay, by, y_rms, len(ym_p)), (c0_a, r0_a, c0_b, r0_b, c1_a, r1_a, c1_b, r1_b)

def read(image, at):
    (ax, bx, x_rms, nx), (ay, by, y_rms, ny), box = _calibrate(image)
    c0, r0, _, _, _, _, c1, r1 = box
    scale = 3
    sub = image[r0:r1, c0:c1]
    H, W, _ = sub.shape
    Hs, Ws = H // scale, W // scale
    small = sub[:Hs * scale, :Ws * scale].reshape(Hs, scale, Ws, scale, 3).mean(axis=(1, 3))
    R, G, B = small[:, :, 0], small[:, :, 1], small[:, :, 2]
    org = ((R > 175) & (G > 90) & (G < 175) & (B < 95)).astype(float)

    kr = 10
    yy, xx = np.ogrid[-kr:kr + 1, -kr:kr + 1]
    sq = (np.maximum(np.abs(xx), np.abs(yy)) <= 8).astype(float)
    sq /= sq.sum()
    sq_sub = np.maximum(np.abs(xx), np.abs(yy)) <= 8

    c_map = signal.fftconvolve(org, sq, mode="same")
    res = org.copy()
    pts = []
    for _ in range(40):
        idx = int(np.argmax(c_map))
        val = float(c_map.flat[idx])
        if val < 0.45:
            break
        ry, rx = divmod(idx, Ws)
        cx = c0 + (rx + 0.5) * scale
        cy = r0 + (ry + 0.5) * scale
        pts.append((cx, cy, val))
        y0, y1 = max(0, ry - kr), min(Hs, ry + kr + 1)
        x0, x1 = max(0, rx - kr), min(Ws, rx + kr + 1)
        ky0, ky1 = y0 - (ry - kr), (2 * kr + 1) - ((ry + kr + 1) - y1)
        kx0, kx1 = x0 - (rx - kr), (2 * kr + 1) - ((rx + kr + 1) - x1)
        removed = np.zeros((y1 - y0, x1 - x0), dtype=float)
        mh = sq_sub[ky0:ky1, kx0:kx1]
        removed[mh] = res[y0:y1, x0:x1][mh]
        res[y0:y1, x0:x1][mh] = 0.0
        delta = signal.convolve2d(removed, sq, mode="full")
        cy0, cx0 = y0 - kr, x0 - kr
        cy1, cx1 = cy0 + delta.shape[0], cx0 + delta.shape[1]
        dy0, dy1 = max(0, -cy0), delta.shape[0] - max(0, cy1 - Hs)
        dx0, dx1 = max(0, -cx0), delta.shape[1] - max(0, cx1 - Ws)
        c_map[max(0, cy0):min(Hs, cy1), max(0, cx0):min(Ws, cx1)] -= delta[dy0:dy1, dx0:dx1]

    pts.sort(key=lambda p: (p[0], -p[1]))
    out = []
    for cx, cy, val in pts:
        x_val = 10.0 ** (ax + bx * cx)
        y_val = 10.0 ** (ay + by * cy)
        half_px = math.sqrt(6.0 ** 2 + (2.0 * y_rms) ** 2)
        y_hi = 10.0 ** (ay + by * (cy - half_px))
        y_lo = 10.0 ** (ay + by * (cy + half_px))
        out.append({
            "x": float(x_val),
            "y": float(y_val),
            "lo": float(y_lo),
            "hi": float(y_hi),
            "px": [float(cx), float(cy)],
            "px_lo": float(cy - half_px),
            "px_hi": float(cy + half_px),
            "note": f"orange filled square (cov={val:.2f}); x_rms={x_rms:.2f}px, y_rms={y_rms:.2f}px",
        })
    return out
