"""Calibrated 1200 dpi reading of Fig. 3D (40.7% development) predicted dashed model profile."""

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


def _calibrate(img):
    gray = img.astype(float).mean(axis=2)
    dark = gray < 80
    H, W = dark.shape
    hr = np.where(dark.mean(axis=1) > 0.5)[0]
    r_top = (int(hr[hr < H // 2][0]), int(hr[hr < H // 2][-1]))
    r_bot = (int(hr[hr > H // 2][0]), int(hr[hr > H // 2][-1]))
    vc = np.where(dark[r_top[0]:r_bot[1] + 1, :].mean(axis=0) > 0.5)[0]
    mid_c = (int(vc[0]) + int(vc[-1])) // 2
    c_left = (int(vc[vc < mid_c][0]), int(vc[vc < mid_c][-1]))
    c_right = (int(vc[vc > mid_c][0]), int(vc[vc > mid_c][-1]))

    xb = _find_ticks(dark[r_bot[0] - 22:r_bot[0] - 6, c_left[1] + 10:c_right[0] - 10], c_left[1] + 10, axis=0)
    yl = _find_ticks(dark[r_top[1] + 25:r_bot[0] - 180, c_left[1] + 4:c_left[1] + 14], r_top[1] + 25, axis=1)
    assert len(xb) == 5, f"expected 5 x-ticks, got {xb}"
    assert len(yl) == 4, f"expected 4 upper y-ticks, got {yl}"

    xa = D.Axis("linear", list(zip(xb, [-0.8, -0.4, 0.0, 0.4, 0.8])))
    ya = D.Axis("linear", list(zip(yl, [100.0, 75.0, 50.0, 25.0])))
    return xa, ya, (c_left, c_right, r_top, r_bot)


def read(image, at):
    xa, ya, (c_left, c_right, r_top, r_bot) = _calibrate(image)
    r0, r1 = r_top[1] + 35, r_bot[0] - 35
    c0, c1 = c_left[1] + 10, c_right[0] - 10
    interior = image[r0:r1, c0:c1]
    sat = interior.max(axis=2).astype(int) - interior.min(axis=2).astype(int)
    col_mask = sat > 55
    lab, n = ndimage.label(col_mask)

    out = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 400:
            continue
        r_c = float(rr.mean()) + 0.5 + r0
        c_c = float(cc.mean()) + 0.5 + c0
        x_val = xa.to_value(c_c)
        if abs(x_val) > 0.72:
            continue
        y_val = ya.to_value(r_c)
        half_px = max(3.0, 0.1 * float(rr.max() - rr.min() + 1)) + ya.residual_px
        a, b = ya.to_value(r_c - half_px), ya.to_value(r_c + half_px)
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": min(a, b),
            "hi": max(a, b),
            "px": [c_c, r_c],
            "px_lo": r_c - half_px,
            "px_hi": r_c + half_px,
            "note": f"magenta model dash centroid ({len(rr)} px; res_x={xa.residual_px:.2f}, res_y={ya.residual_px:.2f} px)",
        })
    return sorted(out, key=lambda d: d["x"])
