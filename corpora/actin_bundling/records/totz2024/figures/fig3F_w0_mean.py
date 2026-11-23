"""Calibrated 1200 dpi reading of Fig. 3F mean (white dot) scaled protoridge amplitudes w_0."""

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
            if 5 <= j - i + 1 <= 25:
                out.append((i + j) / 2.0 + 0.5 + offset)
            i = j + 1
        else:
            i += 1
    return out


def read(image, at):
    gray = image.astype(float).mean(axis=2)
    dark = gray < 80
    H, W = dark.shape
    hr = np.where(dark.mean(axis=1) > 0.5)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(dark[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.5)[0]
    mid_c = (int(vc[0]) + int(vc[-1])) // 2
    c_left = (int(vc[vc < mid_c][0]), int(vc[vc < mid_c][-1]))
    c_right = (int(vc[vc > mid_c][0]), int(vc[vc > mid_c][-1]))

    xt = _find_ticks(dark[r_top[1] + 10:r_top[1] + 30, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yl_in = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_left[1] + 8:c_left[1] + 26], r_top[1] + 10, axis=1)
    assert len(xt) == 5, f"expected 5 x-ticks, got {xt}"
    assert len(yl_in) == 3, f"expected 3 interior y-ticks (150, 100, 50), got {yl_in}"

    r_t_mid = 0.5 * (r_top[0] + r_top[1]) + 0.5
    r_b_mid = 0.5 * (r_bot[0] + r_bot[1]) + 0.5
    ya = D.Axis("linear", list(zip([r_t_mid] + yl_in + [r_b_mid], [200.0, 150.0, 100.0, 50.0, 0.0])))

    times = [float(t) for t in at] if at else [39.6, 40.0, 40.3, 40.7, 41.0]
    err_px = 2.0 + ya.residual_px
    out = []
    for t_val, xc in zip(times, xt):
        ic = int(round(xc))
        strip = image[r_top[1] + 15:r_bot[0] - 15, ic - 30:ic + 31]
        white = (strip[:, :, 0] > 235) & (strip[:, :, 1] > 235) & (strip[:, :, 2] > 235)
        lab, n = ndimage.label(white)
        dot_r, dot_c = None, None
        for k in range(1, n + 1):
            rr, cc = np.nonzero(lab == k)
            if 200 <= len(rr) <= 2500:
                hb = rr.max() - rr.min() + 1
                wb = cc.max() - cc.min() + 1
                if 0.7 <= hb / wb <= 1.4 and abs(cc.mean() - 30) < 10:
                    dot_r = float(rr.mean()) + 0.5 + r_top[1] + 15
                    dot_c = float(cc.mean()) + 0.5 + ic - 30
                    break
        assert dot_r is not None, f"white mean dot not found at t={t_val}%"
        y_val = ya.to_value(dot_r)
        ylo, yhi = sorted((ya.to_value(dot_r - err_px), ya.to_value(dot_r + err_px)))
        out.append({
            "x": t_val,
            "y": y_val,
            "lo": ylo,
            "hi": yhi,
            "px": [dot_c, dot_r],
            "px_lo": dot_r - err_px,
            "px_hi": dot_r + err_px,
            "note": f"white circle mean marker (res_y={ya.residual_px:.2f} px)",
        })
    return out
