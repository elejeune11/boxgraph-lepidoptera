"""Calibrated 1200 dpi reading of Fig. 4C gray 1 nm isocontour L(g)."""

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

    gs = [float(v) for v in at] if at else [0.25, 0.35, 0.45, 0.55, 0.65]
    c0, c1 = int(round(xa.to_pixel(0.355))), int(round(xa.to_pixel(0.455)))
    err_px = 3.0 + xa.residual_px
    out = []
    for g_val in gs:
        r = int(round(ya.to_pixel(g_val)))
        row = image[r, c0:c1]
        g_row = row.mean(axis=1)
        sat_row = row.max(axis=1).astype(float) - row.min(axis=1).astype(float)
        mask = (g_row < 182) & (g_row > 115) & (sat_row < 30)
        s, e = _longest_run(mask, c0)
        assert e - s + 1 >= 10, f"1 nm gray contour not found at g={g_val}%"
        c_mid = 0.5 * (s + e) + 0.5
        L_val = xa.to_value(c_mid)
        lo, hi = sorted((xa.to_value(c_mid - err_px), xa.to_value(c_mid + err_px)))
        out.append({
            "x": g_val,
            "y": L_val,
            "lo": lo,
            "hi": hi,
            "px": [c_mid, float(r) + 0.5],
            "note": f"1 nm gray contour L(g) (run {e-s+1} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return out
