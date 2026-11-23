"""Digitize the 2 red Control markers (minor vs major length, verifying swapped attP40/attP2 labels) in Supplementary Fig. 3c (PDF p. 18, 900 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 8) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read(image, at):
    img = image.astype(float)
    H, W, _ = img.shape
    g = img.mean(axis=2)

    is_ng = (g < 242.0) & (g > 150.0) & (np.abs(img[:, :, 0] - img[:, :, 2]) < 12.0)
    vc = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_ng.sum(axis=0) > H * 0.45)[0], gap=10)]
    hr = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_ng.sum(axis=1) > W * 0.45)[0], gap=10)]
    x_labels = [100.0, 150.0, 200.0, 250.0]
    y_labels = [400.0, 350.0, 300.0, 250.0]
    if len(vc) != 4 or len(hr) != 4:
        raise RuntimeError(f"Expected 4x4 grid in Supp. Fig. 3c, found {len(vc)}x{len(hr)}")

    bx, ax = np.polyfit(vc, x_labels, 1)
    by, ay = np.polyfit(hr, y_labels, 1)
    res_x = float(np.max(np.abs(np.array(vc) - (np.array(x_labels) - ax) / bx)))
    res_y = float(np.max(np.abs(np.array(hr) - (np.array(y_labels) - ay) / by)))
    if res_x > 2.5 or res_y > 2.5:
        raise RuntimeError(f"Grid calibration residual too high: rx={res_x:.2f}px, ry={res_y:.2f}px")

    r0, r1 = int(hr[0]) + 2, int(hr[-1]) - 2
    c0, c1 = int(vc[0]) + 2, int(vc[-1]) - 2
    sub = img[r0:r1, c0:c1]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]

    red = (R - G > 55.0) & (R - B > 55.0) & (R > 180.0)
    opened = ndimage.binary_opening(red, structure=np.ones((15, 15), dtype=bool))
    lab, n = ndimage.label(opened)

    pts = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) < 400:
            continue
        c_abs = float(cc.mean()) + 0.5 + c0
        r_abs = float(rr.mean()) + 0.5 + r0
        x_val = ax + bx * c_abs
        y_val = ay + by * r_abs
        err_px = 2.0 + res_y
        err_y = abs(by) * err_px
        pts.append({
            "x": x_val,
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "px": [c_abs, r_abs],
            "px_lo": r_abs - err_px,
            "px_hi": r_abs + err_px,
            "note": f"red control circle (area={len(rr)}px, rx={res_x:.2f}px, ry={res_y:.2f}px)",
        })
    pts.sort(key=lambda p: p["x"])
    if len(pts) != 2:
        raise RuntimeError(f"Expected 2 red control markers in Supp. Fig. 3c, found {len(pts)}")
    return pts
