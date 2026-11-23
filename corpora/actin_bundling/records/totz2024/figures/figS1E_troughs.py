"""Calibrated 900 dpi reading of Fig. S1E green circle ridge-trough markers."""

import numpy as np
from scipy import ndimage
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis):
    b = strip.mean(axis=axis) > 0.4
    out, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j + 1 < len(b) and b[j + 1]:
                j += 1
            if 1 <= j - i + 1 <= 15:
                out.append((i + j) / 2.0 + 0.5 + offset)
            i = j + 1
        else:
            i += 1
    return out


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 120
    H, W = dark.shape
    hr = np.where(dark.mean(axis=1) > 0.5)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(dark[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.5)[0]
    vgroups = []
    for c in vc:
        if not vgroups or c - vgroups[-1][-1] > 5:
            vgroups.append([int(c)])
        else:
            vgroups[-1].append(int(c))
    c_left = (vgroups[0][0], vgroups[0][-1])
    c_right = (vgroups[-1][0], vgroups[-1][-1])

    xt = _find_ticks(dark[r_top[1] + 10:r_top[1] + 40, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yl = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_left[1] + 10:c_left[1] + 40], r_top[1] + 10, axis=1)
    assert len(xt) == 5 and len(yl) == 3

    c_l_mid = 0.5 * (c_left[0] + c_left[1]) + 0.5
    r_t_mid = 0.5 * (r_top[0] + r_top[1]) + 0.5
    r_b_mid = 0.5 * (r_bot[0] + r_bot[1]) + 0.5
    xa = D.Axis("linear", [(c_l_mid, 0.0)] + list(zip(xt, [2.0, 4.0, 6.0, 8.0, 10.0])))
    ya = D.Axis("linear", list(zip([r_t_mid] + yl + [r_b_mid], [0.4, 0.3, 0.2, 0.1, 0.0])))

    r0, c0 = r_top[1] + 5, c_left[1] + 5
    box = image[r0:r_bot[0] - 5, c0:c_right[0] - 5]
    R, G, B = box[:, :, 0].astype(int), box[:, :, 1].astype(int), box[:, :, 2].astype(int)
    green = (G > 110) & (R < 85) & (B < 85)
    lab, n = ndimage.label(ndimage.binary_dilation(green, iterations=2))

    err_px = 2.0 + ya.residual_px
    out = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 200:
            continue
        r_c = float(rr.mean()) + 0.5 + r0
        c_c = float(cc.mean()) + 0.5 + c0
        x_val = xa.to_value(c_c)
        y_val = ya.to_value(r_c)
        ylo, yhi = sorted((ya.to_value(r_c - err_px), ya.to_value(r_c + err_px)))
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": ylo,
            "hi": yhi,
            "px": [c_c, r_c],
            "px_lo": r_c - err_px,
            "px_hi": r_c + err_px,
            "note": f"green circle ridge trough ({len(rr)} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return sorted(out, key=lambda d: d["x"])
