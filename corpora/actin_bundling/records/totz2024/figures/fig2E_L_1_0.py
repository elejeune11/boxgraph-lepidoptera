"""Calibrated 1200 dpi reading of Fig. 2E (Bundle distance 1.0 L, medium purple curve)."""

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
    hr = np.where(dark.mean(axis=1) > 0.5)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(dark[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.5)[0]
    mid_c = (int(vc[0]) + int(vc[-1])) // 2
    c_left = (int(vc[vc < mid_c][0]), int(vc[vc < mid_c][-1]))
    c_right = (int(vc[vc > mid_c][0]), int(vc[vc > mid_c][-1]))

    xt = _find_ticks(dark[r_top[1] + 5:r_top[1] + 16, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yl = _find_ticks(dark[r_top[1] + 10:r_bot[0] - 10, c_left[1] + 5:c_left[1] + 16], r_top[1] + 10, axis=1)
    assert len(xt) == 5 and len(yl) == 4
    r_b_mid = 0.5 * (r_bot[0] + r_bot[1]) + 0.5
    xa = D.Axis("linear", list(zip(xt, [-0.5, -0.25, 0.0, 0.25, 0.5])))
    ya = D.Axis("linear", list(zip(yl + [r_b_mid], [400.0, 300.0, 200.0, 100.0, 0.0])))

    r_min = int(round(ya.to_pixel(250.0)))
    r_max = r_bot[0] - 4
    xs = [float(v) for v in at] if at else [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]
    err_px = 2.0 + ya.residual_px
    out = []
    for x_val in xs:
        c = int(round(xa.to_pixel(x_val)))
        col = image[r_min:r_max, c]
        R, G, B = col[:, 0].astype(int), col[:, 1].astype(int), col[:, 2].astype(int)
        mask = (R > 95) & (R < 175) & (G < 50) & (B > 95) & (B < 175)
        rr = np.nonzero(mask)[0] + r_min
        assert len(rr) >= 5, f"medium purple curve not found at x/L={x_val}"
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
            "note": f"medium purple 1.0L curve (stroke {len(rr)} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return out
