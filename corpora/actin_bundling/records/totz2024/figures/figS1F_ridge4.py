"""Calibrated 1200 dpi reading of Fig. S1F (ridge 4) segmented & rotated red mean height profile."""

import numpy as np
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis):
    b = strip.mean(axis=axis) > 0.35
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
    dark = gray < 130
    H, W = dark.shape
    vc = np.where(dark.sum(axis=0) > 0.5 * H)[0]
    vgroups = []
    for c in vc:
        if not vgroups or c - vgroups[-1][-1] > 10:
            vgroups.append([int(c)])
        else:
            vgroups[-1].append(int(c))
    c_left = (vgroups[0][0], vgroups[0][-1])
    c_right = (vgroups[-1][0], vgroups[-1][-1])
    hr = np.where(dark[:, c_left[0]:c_right[1] + 1].mean(axis=1) > 0.6)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))

    xt = _find_ticks(dark[r_top[1] + 3:r_top[1] + 12, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yr = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_right[0] - 12:c_right[0] - 3], r_top[1] + 10, axis=1)
    assert len(xt) == 3 and len(yr) == 3, f"ticks xt={xt}, yr={yr}"

    clm = 0.5 * (c_left[0] + c_left[1]) + 0.5
    crm = 0.5 * (c_right[0] + c_right[1]) + 0.5
    rtm = 0.5 * (r_top[0] + r_top[1]) + 0.5
    xa = D.Axis("linear", [(clm, 0.0)] + list(zip(xt, [0.5, 1.0, 1.5])) + [(crm, 2.0)])
    ya = D.Axis("linear", list(zip([rtm] + yr, [0.15, 0.10, 0.05, 0.0])))

    r0, r1 = r_top[1] + 5, r_bot[0] - 5
    xs = [float(v) for v in at] if at else [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    err_px = 2.5 + ya.residual_px
    out = []
    for x_val in xs:
        c = int(round(xa.to_pixel(x_val)))
        col = image[r0:r1, c]
        R, G, B = col[:, 0].astype(int), col[:, 1].astype(int), col[:, 2].astype(int)
        rr = np.nonzero((R > 165) & (G < 85) & (B < 85))[0] + r0
        assert len(rr) >= 3, f"red curve not found at x={x_val}"
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
            "note": f"red rotated mean ridge 4 (stroke {len(rr)} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return out
