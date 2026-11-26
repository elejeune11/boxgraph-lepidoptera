"""Calibrated curve reader for fig3d_asymptote_decoupled (Fig. 3, p. 19)."""

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

COLOR_KIND = 'dash'


def _guide(n):
    return (0.05 * n, 0.095 * n)

def read(image, at):
    (ax, bx, x_rms, nx), (ay, by, y_rms, ny), box = _calibrate(image)
    _, _, c0, r0, c1, r1, _, _ = box
    img_f = image.astype(float)
    R_all, G_all, B_all = img_f[:, :, 0], img_f[:, :, 1], img_f[:, :, 2]

    # Mask solid marker halos so anti-aliased marker rims never pull thin curves off-center
    solid_red = (R_all > 135) & (G_all < 85) & (B_all < 85) & ((R_all - G_all) > 75)
    solid_org = (R_all > 175) & (G_all > 90) & (G_all < 175) & (B_all < 95)
    marker_halo = ndimage.binary_dilation(solid_red | solid_org, iterations=18)

    out = []
    for x_val in at:
        n_val = float(x_val)
        col = int(round((math.log10(n_val) - ax) / bx))
        y_lo_g, y_hi_g = _guide(n_val)
        ra = max(r0 + 18, int(round((math.log10(y_hi_g) - ay) / by)))
        rb = min(r1 - 18, int(round((math.log10(y_lo_g) - ay) / by)))
        if rb <= ra + 4:
            out.append({"x": n_val, "y": None, "note": "guide band outside plot box"})
            continue

        if COLOR_KIND == "dash":
            ca = max(c0 + 10, col - 35)
            cb = min(c1 - 10, col + 36)
            patch = img_f[ra:rb, ca:cb]
            pgray = patch.mean(axis=2)
            psat = patch.max(axis=2) - patch.min(axis=2)
            pmask = (pgray < 215) & (psat < 18) & (~marker_halo[ra:rb, ca:cb])
            rr, cc = np.where(pmask)
            if len(rr) < 12:
                out.append({"x": n_val, "y": None, "note": "dashed line not found"})
                continue
            c_rel = cc + ca - col
            m, b_lin = np.polyfit(c_rel, rr, 1)
            for _ in range(2):
                res = np.abs(rr - (b_lin + m * c_rel))
                keep = res < 12.0
                if keep.sum() >= 8:
                    m, b_lin = np.polyfit(c_rel[keep], rr[keep], 1)
            r_c = ra + float(b_lin) + 0.5
            slope_err = abs(float(m)) * (2.0 * x_rms)
            half_px = math.sqrt(7.0 ** 2 + (2.0 * y_rms) ** 2 + slope_err ** 2)
        else:
            strip = img_f[ra:rb, col - 2:col + 3, :].mean(axis=1)
            halo_col = marker_halo[ra:rb, col - 2:col + 3].any(axis=1)
            R, G, B = strip[:, 0], strip[:, 1], strip[:, 2]
            gray = strip.mean(axis=1)
            if COLOR_KIND == "red":
                score = np.where(
                    (R - G > 10) & (R - B > 4) & (G > 90) & (gray < 235) & (~halo_col),
                    (R - G) + (245.0 - G) * 0.3,
                    0.0,
                )
            elif COLOR_KIND == "blue":
                score = np.where(
                    (B - G >= -1.5) & (B - R >= -3.0) & (gray < 215) & (gray > 75) & ((R - G) < 10) & (~halo_col),
                    (240.0 - gray),
                    0.0,
                )
            elif COLOR_KIND == "green":
                score = np.where(
                    (G - R > 6) & (G - B > 4) & (gray < 235) & (~halo_col),
                    (G - R) + (245.0 - R) * 0.2,
                    0.0,
                )
            elif COLOR_KIND == "orange":
                score = np.where(
                    (R - B > 10) & (G - B > 8) & (G > 170) & (R > 180) & (B < 230) & (~halo_col),
                    (R - B) + (G - B),
                    0.0,
                )
            else:
                raise ValueError(COLOR_KIND)

            idx = np.where(score > 0)[0]
            if len(idx) == 0:
                out.append({"x": n_val, "y": None, "note": f"{COLOR_KIND} curve not found"})
                continue
            gs = [[int(idx[0])]]
            for v in idx[1:]:
                iv = int(v)
                if iv <= gs[-1][-1] + 2:
                    gs[-1].append(iv)
                else:
                    gs.append([iv])
            thin = [g for g in gs if len(g) <= 28]
            if not thin:
                out.append({"x": n_val, "y": None, "note": "only thick feature in guide band"})
                continue
            best_g = max(thin, key=lambda g: float(score[g[0]:g[-1] + 1].sum()))
            w = score[best_g[0]:best_g[-1] + 1]
            r_c = ra + best_g[0] + float(np.sum(np.arange(len(best_g)) * w) / w.sum()) + 0.5
            stroke_half = max(6.0, 0.5 * len(best_g))
            half_px = math.sqrt(stroke_half ** 2 + (2.0 * y_rms) ** 2 + (1.5 * 2.0 * x_rms) ** 2)

        y_val = 10.0 ** (ay + by * r_c)
        y_hi = 10.0 ** (ay + by * (r_c - half_px))
        y_lo = 10.0 ** (ay + by * (r_c + half_px))
        out.append({
            "x": n_val,
            "y": float(y_val),
            "lo": float(y_lo),
            "hi": float(y_hi),
            "px": [float(col), float(r_c)],
            "px_lo": float(r_c - half_px),
            "px_hi": float(r_c + half_px),
            "note": f"x_rms={x_rms:.2f}px ({nx} ticks), y_rms={y_rms:.2f}px ({ny} ticks), half={half_px:.1f}px",
        })
    return out
