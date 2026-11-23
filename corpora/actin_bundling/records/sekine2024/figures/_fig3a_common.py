"""Shared calibration and colored-marker reader for Fig. 3a (RNAi/mutant screen scatter, PDF p. 5, 900 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 8) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read_fig3a_markers(image: np.ndarray, color_group: str):
    """Calibrate Fig. 3a from its 4x5 rectangular gray grid and read red, blue, or green circular markers."""
    img = image.astype(float)
    H, W, _ = img.shape
    gray = img.mean(axis=2)
    is_grid = (gray < 235.0) & (gray > 140.0) & (np.abs(img[:, :, 0] - img[:, :, 2]) < 10.0)

    hr = [float(g.mean()) + 0.5 for g in _split_runs(np.where(is_grid.sum(axis=1) > W * 0.55)[0], gap=10)]
    vc = [float(g.mean()) + 0.5 for g in _split_runs(np.where(is_grid.sum(axis=0) > H * 0.55)[0], gap=10)]
    x_labels = [40.0, 50.0, 60.0, 70.0]
    y_labels = [400.0, 360.0, 320.0, 280.0, 240.0]
    if len(vc) != len(x_labels) or len(hr) != len(y_labels):
        raise RuntimeError(f"Expected 4x5 grid, found {len(vc)}x{len(hr)}")

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

    if color_group == "red":
        mask = (R - G > 55.0) & (R - B > 55.0) & (R > 180.0)
    elif color_group == "blue":
        mask = (B - R > 30.0) & (B > 160.0)
    elif color_group == "green":
        mask = (G - R > 18.0) & (G - B > 25.0) & (G > 140.0)
    else:
        raise ValueError(f"Unknown color_group: {color_group}")

    # Opening by 15x15 eliminates any colored text (e.g. red Control labels) while preserving filled circle interiors
    opened = ndimage.binary_opening(mask, structure=np.ones((15, 15), dtype=bool))
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
            "note": f"{color_group} circle marker (area={len(rr)}px, rx={res_x:.2f}px, ry={res_y:.2f}px)",
        })
    pts.sort(key=lambda p: (round(p["x"], 1), p["y"]))
    return pts
