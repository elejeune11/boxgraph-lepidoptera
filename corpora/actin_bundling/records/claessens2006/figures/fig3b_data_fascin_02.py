"""Calibrated blue-open-circle marker reader for fig3b_data_fascin_02 (Fig. 3, p. 19)."""

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
    pad = 36
    scale = 3
    y_off, x_off = max(0, r0 - pad), max(0, c0 - pad)
    sub = image[y_off:min(image.shape[0], r1 + pad), x_off:min(image.shape[1], c1 + pad)]
    H, W, _ = sub.shape
    Hs, Ws = H // scale, W // scale
    small = sub[:Hs * scale, :Ws * scale].reshape(Hs, scale, Ws, scale, 3).mean(axis=(1, 3))
    R, G, B = small[:, :, 0], small[:, :, 1], small[:, :, 2]
    gray = small.mean(axis=2)

    white = (R > 215) & (G > 215) & (B > 215)
    lbl, _ = ndimage.label(white)
    objs = ndimage.find_objects(lbl)

    blue_ring = ((gray < 195) & ((R - B) < 35) & ((G - R) < 10) & (B >= G - 1)).astype(float)
    dark_navy = ((gray < 130) & ((R - B) < 35) & ((G - R) < 10)).astype(float)
    yy, xx = np.ogrid[-12:13, -12:13]
    rr_k = np.sqrt(xx * xx + yy * yy)
    k_ring = ((rr_k >= 6.8) & (rr_k <= 11.2)).astype(float)
    k_ring /= k_ring.sum()
    ring_map = signal.fftconvolve(blue_ring, k_ring, mode="same")
    dark_map = signal.fftconvolve(dark_navy, k_ring, mode="same")
    k_in = (rr_k <= 5.5).astype(float)
    k_in /= k_in.sum()
    in_map = signal.fftconvolve(white.astype(float), k_in, mode="same")

    pts = []
    for i, sl in enumerate(objs, 1):
        if sl is None:
            continue
        comp = lbl[sl] == i
        area = int(comp.sum())
        h_box = sl[0].stop - sl[0].start
        w_box = sl[1].stop - sl[1].start
        if not (18 <= area <= 260) or w_box > 17 or h_box > 17:
            continue
        ys, xs = np.where(comp)
        ry_i = int(round(sl[0].start + float(ys.mean())))
        rx_i = int(round(sl[1].start + float(xs.mean())))
        y0, y1 = max(0, ry_i - 5), min(Hs, ry_i + 6)
        x0, x1 = max(0, rx_i - 5), min(Ws, rx_i + 6)
        cand = ring_map[y0:y1, x0:x1] * (in_map[y0:y1, x0:x1] >= 0.35)
        if cand.max() < 0.28:
            continue
        dy_m, dx_m = np.unravel_index(int(np.argmax(cand)), cand.shape)
        by_c, bx_c = y0 + dy_m, x0 + dx_m
        if dark_map[by_c, bx_c] < 0.16:
            continue
        # Verify ring is bluish/navy (not green square) using non-red dark ring pixels
        ry0, ry1 = max(0, by_c - 11), min(Hs, by_c + 12)
        rx0, rx1 = max(0, bx_c - 11), min(Ws, bx_c + 12)
        py, px_g = np.ogrid[ry0 - by_c:ry1 - by_c, rx0 - bx_c:rx1 - bx_c]
        d = np.sqrt(px_g * px_g + py * py)
        Rp, Gp, Bp, grp = R[ry0:ry1, rx0:rx1], G[ry0:ry1, rx0:rx1], B[ry0:ry1, rx0:rx1], gray[ry0:ry1, rx0:rx1]
        rmask = (d >= 6.8) & (d <= 11.2) & (grp < 185) & ((Rp - Bp) < 30)
        if rmask.sum() < 18:
            continue
        if float((Gp[rmask] - Rp[rmask]).mean()) > 5.0:
            continue
        if float((Bp[rmask] - Gp[rmask]).mean()) < 1.2:
            continue
        cx = x_off + (bx_c + 0.5) * scale
        cy = y_off + (by_c + 0.5) * scale
        if c0 - 15 <= cx <= c1 + 15 and r0 - 15 <= cy <= r1 + 15:
            if all((cx - p[0]) ** 2 + (cy - p[1]) ** 2 > 18 ** 2 for p in pts):
                pts.append((cx, cy, area, float(ring_map[by_c, bx_c])))

    pts.sort(key=lambda p: (p[0], -p[1]))
    out = []
    for cx, cy, area, rscore in pts:
        x_val = 10.0 ** (ax + bx * cx)
        y_val = 10.0 ** (ay + by * cy)
        loc_err = 6.0 if area >= 100 else 9.0
        half_px = math.sqrt(loc_err ** 2 + (2.0 * y_rms) ** 2)
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
            "note": f"blue open circle (area={area}, ring={rscore:.2f}); x_rms={x_rms:.2f}px, y_rms={y_rms:.2f}px",
        })
    return out
