"""Calibrated 1200 dpi reading of Fig. 1F (40.0% pupal development) 1D in vivo centerline height profile."""

import numpy as np
from boxgraph import digitize as D


def _find_ticks(strip, offset, axis, thresh=0.35):
    b = strip.mean(axis=axis) > thresh
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
    dark = gray < 130
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
    c_left = (vgroups[-2][0], vgroups[-2][-1])
    c_right = (vgroups[-1][0], vgroups[-1][-1])

    xb = _find_ticks(dark[r_bot[0] - 7:r_bot[0] - 2, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yr = _find_ticks(dark[r_top[1] + 5:r_bot[0] - 5, c_right[0] - 8:c_right[0] - 2], r_top[1] + 5, axis=1)
    assert len(xb) == 6, f"expected 6 x-ticks (0..5), got {xb}"
    assert len(yr) == 3, f"expected 3 y-ticks (100, 50, 0), got {yr}"

    xa = D.Axis("linear", list(zip(xb, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])))
    ya = D.Axis("linear", list(zip(yr, [100.0, 50.0, 0.0])))

    r0 = r_top[1] + 8
    r1 = r_bot[0] - 8
    xs = [float(v) for v in at] if at else [0.25 * k for k in range(1, 20)]
    err_px = 2.0 + ya.residual_px
    out = []
    for x_val in xs:
        c = int(round(xa.to_pixel(x_val)))
        col = gray[r0:r1, c]
        rr = np.nonzero(col < 110)[0] + r0
        assert len(rr) >= 3, f"black curve not found at x={x_val}"
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
            "note": f"in vivo centerline height profile (stroke {len(rr)} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return out
