"""Digitize all 10 markers (8 CGMD simulations + St.15, St.16) of major vs minor length in Supplementary Fig. 6b (PDF p. 21, 1200 dpi)."""

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
    vc = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_ng.sum(axis=0) > H * 0.25)[0], gap=10)]
    hr = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_ng.sum(axis=1) > W * 0.20)[0], gap=10)]
    x_labels = [100.0, 150.0, 200.0, 250.0, 300.0]
    y_labels = [450.0, 400.0, 350.0, 300.0, 250.0]
    if len(vc) != 5 or len(hr) != 5:
        raise RuntimeError(f"Expected 5x5 grid in Supp. Fig. 6b, found {len(vc)}x{len(hr)}")

    bx, ax = np.polyfit(vc, x_labels, 1)
    by, ay = np.polyfit(hr, y_labels, 1)
    res_x = float(np.max(np.abs(np.array(vc) - (np.array(x_labels) - ax) / bx)))
    res_y = float(np.max(np.abs(np.array(hr) - (np.array(y_labels) - ay) / by)))
    if res_x > 5.0 or res_y > 5.0:
        raise RuntimeError(f"Grid calibration residual too high: rx={res_x:.2f}px, ry={res_y:.2f}px")

    r0, r1 = int(hr[0]) + 4, int(hr[-1]) - 4
    c0, c1 = int(vc[0]) + 4, int(vc[-1]) - 4
    sub = img[r0:r1, c0:c1]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    gsub = sub.mean(axis=2)

    specs = [
        ("red_4", ((R > 205) & (R - G > 55) & (B > 60) & (G < 185)) | ((R > 200) & (G < 80) & (B < 80)), 15, 700),
        ("orange_2", (R > 205) & (G > 95) & (G < 205) & (B < 145) & (R - B > 95) & (R - G > 25), 11, 500),
        ("yellow_3", (R > 200) & (G > 170) & (B < 165) & (G - B > 50) & (R - G < 55), 15, 700),
        ("blue_1", (B > 170) & (B - R > 35) & (G > 120), 15, 700),
        ("gray_exp", (gsub < 205) & (gsub > 95) & (np.abs(R - B) < 12) & (np.abs(R - G) < 12), 15, 650),
    ]

    pts = []
    for label_kind, mask, ksz, min_area in specs:
        closed = ndimage.binary_closing(mask, structure=np.ones((5, 5), dtype=bool))
        opened = ndimage.binary_opening(closed, structure=np.ones((ksz, ksz), dtype=bool))
        lab, n = ndimage.label(opened)
        for k in range(1, n + 1):
            rr, cc = np.nonzero(lab == k)
            if len(rr) < min_area:
                continue
            h_bb = rr.max() - rr.min() + 1
            w_bb = cc.max() - cc.min() + 1
            if max(h_bb, w_bb) / max(1, min(h_bb, w_bb)) > 1.8:
                continue
            c_abs = 0.5 * (float(cc.min()) + float(cc.max())) + 0.5 + c0
            r_abs = 0.5 * (float(rr.min()) + float(rr.max())) + 0.5 + r0
            if any(abs(p["px"][0] - c_abs) < 18.0 and abs(p["px"][1] - r_abs) < 18.0 for p in pts):
                continue
            x_val = ax + bx * c_abs
            y_val = ay + by * r_abs
            err_px = 2.5 + res_y
            err_y = abs(by) * err_px
            pts.append({
                "x": x_val,
                "y": y_val,
                "lo": y_val - err_y,
                "hi": y_val + err_y,
                "px": [c_abs, r_abs],
                "px_lo": r_abs - err_px,
                "px_hi": r_abs + err_px,
                "note": f"{label_kind} marker (area={len(rr)}px, rx={res_x:.2f}px, ry={res_y:.2f}px)",
            })

    pts.sort(key=lambda p: (round(p["x"], 1), p["y"]))
    if len(pts) != 10:
        raise RuntimeError(f"Expected 10 markers in Supp. Fig. 6b, found {len(pts)}: {[(round(p['x'],1), round(p['y'],1), p['note']) for p in pts]}")
    return pts
