"""Calibrated 1200 dpi reading of Fig. 4C gray 10 nm isocontour g(L)."""

import numpy as np
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis):
    b = strip.mean(axis=axis) > 0.6
    out, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j + 1 < len(b) and b[j + 1]:
                j += 1
            if 4 <= j - i + 1 <= 25:
                out.append((i + j) / 2.0 + 0.5 + offset)
            i = j + 1
        else:
            i += 1
    return out


def _longest_run(mask, offset):
    best_s, best_e = 0, -1
    i = 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j + 1 < len(mask) and mask[j + 1]:
                j += 1
            if j - i > best_e - best_s:
                best_s, best_e = i, j
            i = j + 1
        else:
            i += 1
    return best_s + offset, best_e + offset


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 80
    H, W = dark.shape
    sub = dark[:, :int(W * 0.87)]
    hr = np.where(sub.mean(axis=1) > 0.6)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(sub[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.6)[0]
    mid_c = (int(vc[0]) + int(vc[-1])) // 2
    c_left = (int(vc[vc < mid_c][0]), int(vc[vc < mid_c][-1]))
    c_right = (int(vc[vc > mid_c][0]), int(vc[vc > mid_c][-1]))

    xb_long = _find_ticks(dark[r_bot[0] - 42:r_bot[0] - 36, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yl = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_left[1] + 6:c_left[1] + 26], r_top[1] + 10, axis=1)
    assert len(xb_long) == 3 and len(yl) == 3
    c_r_mid = 0.5 * (c_right[0] + c_right[1]) + 0.5
    r_b_mid = 0.5 * (r_bot[0] + r_bot[1]) + 0.5
    xa = D.Axis("linear", [(xb_long[0], 0.5), (xb_long[1], 1.0), (xb_long[2], 1.5), (c_r_mid, 2.0)])
    ya = D.Axis("linear", [(yl[0], 0.6), (yl[1], 0.4), (yl[2], 0.2), (r_b_mid, 0.0)])

    Ls = [float(v) for v in at] if at else [0.51, 0.62, 0.66, 0.70, 0.75]
    r0, r1 = r_top[1] + 15, r_bot[0] - 15
    err_px = 4.0 + ya.residual_px
    out = []
    for L_val in Ls:
        c = int(round(xa.to_pixel(L_val)))
        col = image[r0:r1, c]
        g_col = col.mean(axis=1)
        sat_col = col.max(axis=1).astype(float) - col.min(axis=1).astype(float)
        mask = (g_col < 182) & (g_col > 105) & (sat_col < 30)
        s, e = _longest_run(mask, r0)
        assert e - s + 1 >= 15, f"10 nm gray contour not found at L={L_val}"
        r_mid = 0.5 * (s + e) + 0.5
        g_val = ya.to_value(r_mid)
        glo, ghi = sorted((ya.to_value(r_mid - err_px), ya.to_value(r_mid + err_px)))
        out.append({
            "x": L_val,
            "y": g_val,
            "lo": glo,
            "hi": ghi,
            "px": [float(c) + 0.5, r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"10 nm gray contour (run {e-s+1} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return out
