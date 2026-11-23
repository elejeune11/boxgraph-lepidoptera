"""Calibrated 1200 dpi reading of Fig. S1J rescaled all-ridge mean (red) and +/-1 SD (blue) height profile."""

import numpy as np
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis):
    b = strip.mean(axis=axis) > 0.5
    out, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j + 1 < len(b) and b[j + 1]:
                j += 1
            if 2 <= j - i + 1 <= 15:
                out.append((i + j) / 2.0 + 0.5 + offset)
            i = j + 1
        else:
            i += 1
    return out


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100
    sub = dark[200:, :]
    hr = np.where(sub.mean(axis=1) > 0.4)[0] + 200
    r_top = (int(hr[hr < 700][0]), int(hr[hr < 700][-1]))
    r_bot = (int(hr[hr > 1000][0]), int(hr[hr > 1000][-1]))
    vc = np.where(dark[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.5)[0]
    vgroups = []
    for c in vc:
        if not vgroups or c - vgroups[-1][-1] > 5:
            vgroups.append([int(c)])
        else:
            vgroups[-1].append(int(c))
    c_left = (vgroups[0][0], vgroups[0][-1])
    c_right = (vgroups[-1][0], vgroups[-1][-1])

    xt = _find_ticks(dark[r_top[1] + 3:r_top[1] + 15, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yl = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_left[1] + 3:c_left[1] + 15], r_top[1] + 10, axis=1)
    assert len(xt) == 1 and len(yl) >= 3
    r_t_mid = 0.5 * (r_top[0] + r_top[1]) + 0.5
    c_l_mid = 0.5 * (c_left[0] + c_left[1]) + 0.5
    c_r_mid = 0.5 * (c_right[0] + c_right[1]) + 0.5
    xa = D.Axis("linear", [(c_l_mid, 0.0), (xt[0], 0.8), (c_r_mid, 1.6)])
    # Top spine is 0.20 um; upper 3 interior left ticks are 0.15, 0.10, 0.05 um (0.00 tick is occluded by curves at left edge)
    ya = D.Axis("linear", list(zip([r_t_mid] + yl[:3], [0.20, 0.15, 0.10, 0.05])))

    xs = [float(v) for v in at] if at else [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
    r0, r1 = r_top[1] + 10, r_bot[0] - 10
    err_px = 2.5 + ya.residual_px
    out = []
    for x_val in xs:
        c = int(round(xa.to_pixel(x_val)))
        col = image[r0:r1, c]
        R, G, B = col[:, 0].astype(int), col[:, 1].astype(int), col[:, 2].astype(int)
        is_red = (R > 160) & (G < 80) & (B < 80)
        is_blue = (B > 140) & (R < 80) & (G > 90)
        rr_red = np.nonzero(is_red)[0] + r0
        rr_blue = np.nonzero(is_blue)[0] + r0
        assert len(rr_red) >= 5, f"red mean curve not found at x={x_val}"
        r_mean = float(rr_red.mean()) + 0.5
        up = rr_blue[rr_blue < r_mean - 25]
        lo = rr_blue[rr_blue > r_mean + 25]
        assert len(up) >= 5 and len(lo) >= 5, f"blue SD curves not found at x={x_val}"
        r_up = float(up.mean()) + 0.5
        r_lo = float(lo.mean()) + 0.5

        y_val = ya.to_value(r_mean)
        ylo, yhi = sorted((ya.to_value(r_mean - err_px), ya.to_value(r_mean + err_px)))
        bar_lo, bar_hi = sorted((ya.to_value(r_lo), ya.to_value(r_up)))
        bar_lo_out = min(ya.to_value(r_lo + err_px), bar_lo)
        bar_hi_out = max(ya.to_value(r_up - err_px), bar_hi)
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": ylo,
            "hi": yhi,
            "bar_lo": bar_lo,
            "bar_hi": bar_hi,
            "bar_lo_outer": bar_lo_out,
            "bar_hi_outer": bar_hi_out,
            "px": [float(c) + 0.5, r_mean],
            "px_lo": r_up,
            "px_hi": r_lo,
            "note": f"red mean & blue +/-1 SD curves (res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return out
