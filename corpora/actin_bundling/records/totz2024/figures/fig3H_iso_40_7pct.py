"""Calibrated 900 dpi reading of Fig. 3H pink/magenta (40.7%) isocontour g(Y)."""

import numpy as np
from scipy import ndimage
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis):
    b = strip.mean(axis=axis) > 0.5
    out, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j + 1 < len(b) and b[j + 1]:
                j += 1
            if 2 <= j - i + 1 <= 20:
                out.append((i + j) / 2.0 + 0.5 + offset)
            i = j + 1
        else:
            i += 1
    return out


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 80
    H, W = dark.shape
    sub = dark[:, :int(W * 0.85)]
    hr = np.where(sub.mean(axis=1) > 0.6)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(sub[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.6)[0]
    mid_c = (int(vc[0]) + int(vc[-1])) // 2
    c_left = (int(vc[vc < mid_c][0]), int(vc[vc < mid_c][-1]))
    c_right = (int(vc[vc > mid_c][0]), int(vc[vc > mid_c][-1]))

    xt_long = _find_ticks(dark[r_top[1] + 18:r_top[1] + 28, c_left[1] + 5:c_right[0] - 5], c_left[1] + 5, axis=0)
    yl = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_left[1] + 10:c_left[1] + 28], r_top[1] + 10, axis=1)
    assert len(xt_long) == 2 and len(yl) == 3

    c_l_mid = 0.5 * (c_left[0] + c_left[1]) + 0.5
    c_r_mid = 0.5 * (c_right[0] + c_right[1]) + 0.5
    r_b_mid = 0.5 * (r_bot[0] + r_bot[1]) + 0.5
    xa = D.Axis("log10", [(c_l_mid, 1e5), (xt_long[0], 1e6), (xt_long[1], 1e7), (c_r_mid, 1e8)])
    ya = D.Axis("linear", [(yl[0], 0.6), (yl[1], 0.4), (yl[2], 0.2), (r_b_mid, 0.0)])

    r0, c0 = r_top[1] + 5, c_left[1] + 5
    box_img = image[r0:r_bot[0] - 5, c0:c_right[0] - 5]
    R, G, B = box_img[:, :, 0].astype(int), box_img[:, :, 1].astype(int), box_img[:, :, 2].astype(int)
    mask = (R > 195) & (G > 75) & (G < 130) & (B > 135) & (B < 190)
    lab, n = ndimage.label(mask)

    err_px = 3.0 + ya.residual_px
    out = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 150:
            continue
        r_c = float(rr.mean()) + 0.5 + r0
        c_c = float(cc.mean()) + 0.5 + c0
        Y_val = xa.to_value(c_c)
        g_val = ya.to_value(r_c)
        glo, ghi = sorted((ya.to_value(r_c - err_px), ya.to_value(r_c + err_px)))
        out.append({
            "x": Y_val,
            "y": g_val,
            "lo": glo,
            "hi": ghi,
            "px": [c_c, r_c],
            "px_lo": r_c - err_px,
            "px_hi": r_c + err_px,
            "note": f"magenta 40.7% dash centroid ({len(rr)} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return sorted(out, key=lambda d: d["x"])
