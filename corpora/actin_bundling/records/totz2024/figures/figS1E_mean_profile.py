"""Calibrated 900 dpi reading of Fig. S1E multi-cross-sectional red mean height profile."""

import numpy as np
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

    r0, r1 = r_top[1] + 10, r_bot[0] - 10
    xs = [float(v) for v in at] if at else [0.5 * k for k in range(1, 20)]
    err_px = 2.5 + ya.residual_px
    out = []
    for x_val in xs:
        c = int(round(xa.to_pixel(x_val)))
        # Probe a small 5-px column neighborhood in case a blue/green ring edge touches c
        patch = image[r0:r1, c - 2:c + 3]
        R, G, B = patch[:, :, 0].astype(int), patch[:, :, 1].astype(int), patch[:, :, 2].astype(int)
        is_red = (R > 165) & (G < 85) & (B < 85)
        rr = np.nonzero(is_red.any(axis=1))[0] + r0
        assert len(rr) >= 3, f"red mean curve not found at x={x_val}"
        r_mid = float(rr.mean()) + 0.5
        y_val = ya.to_value(r_mid)
        ylo, yhi = sorted((ya.to_value(r_mid - err_px), ya.to_value(r_mid + err_px)))
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": ylo,
            "hi": yhi,
            "px": [float(c) + 0.5, r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"red mean ROI cross-section (stroke {len(rr)} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return out
