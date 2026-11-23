"""Digitize all 27 blue bars of Average major length (nm) in Supplementary Fig. 3b (PDF p. 18, 900 dpi)."""

from __future__ import annotations

import numpy as np


def _split_runs(idx: np.ndarray, gap: int = 6) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read(image, at):
    img = image.astype(float)
    H, W, _ = img.shape
    g = img.mean(axis=2)

    # 10 horizontal dotted gridlines (450, 400, ..., 0 nm)
    is_ng = (g < 245.0) & (g > 160.0) & (np.abs(img[:, :, 0] - img[:, :, 2]) < 10.0)
    hr = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_ng.sum(axis=1) > W * 0.35)[0], gap=8)]
    y_labels = [450.0, 400.0, 350.0, 300.0, 250.0, 200.0, 150.0, 100.0, 50.0, 0.0]
    if len(hr) != len(y_labels):
        raise RuntimeError(f"Expected 10 horizontal gridlines in Supp. Fig. 3b, found {len(hr)}: {hr}")

    by, ay = np.polyfit(hr, y_labels, 1)
    res_y = float(np.max(np.abs(np.array(hr) - (np.array(y_labels) - ay) / by)))
    if res_y > 3.0:
        raise RuntimeError(f"Y calibration residual too high: ry={res_y:.2f}px")

    r_top_box = int(hr[0])
    r_bot_box = int(hr[-1])
    sub = img[r_top_box:r_bot_box, :]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    is_blue = (B - R > 45.0) & (B > 130.0)

    # Each bar is ~38 px wide and > 400 px tall
    bar_groups = [gr for gr in _split_runs(np.where(is_blue.sum(axis=0) > 250)[0], gap=8) if len(gr) >= 18]
    if len(bar_groups) != 27:
        raise RuntimeError(f"Expected 27 blue bars in Supp. Fig. 3b, found {len(bar_groups)}")

    out = []
    for i, gr in enumerate(bar_groups):
        c_mid = 0.5 * (float(gr[0]) + float(gr[-1])) + 0.5
        # Sample interior columns of bar (avoiding left/right 3px anti-aliased edges)
        c0 = int(gr[0]) + 3
        c1 = int(gr[-1]) - 3
        tops = []
        for cc in range(c0, c1 + 1):
            rr = np.nonzero(is_blue[:, cc])[0]
            if len(rr):
                tops.append(float(rr.min()) + r_top_box)
        r_bar_top = float(np.median(tops))
        y_val = ay + by * r_bar_top
        err_px = 1.5 + res_y
        err_y = abs(by) * err_px
        out.append({
            "x": float(i + 1),
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "px": [c_mid, r_bar_top],
            "px_lo": r_bar_top - err_px,
            "px_hi": r_bar_top + err_px,
            "note": f"bar {i+1} (w={len(gr)}px, ry={res_y:.2f}px)",
        })
    return out
