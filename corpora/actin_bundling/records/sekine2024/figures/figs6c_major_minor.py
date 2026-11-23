"""Digitize all 7 markers (1-M, 1-H circles + a3-1, a4-1, a5-2, a6-1, a6-3 X's) in Supplementary Fig. 6c (PDF p. 21, 1200 dpi)."""

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
    x_labels = [150.0, 200.0, 250.0, 300.0]
    y_labels = [400.0, 350.0, 300.0, 250.0]
    if len(vc) != 4 or len(hr) != 4:
        raise RuntimeError(f"Expected 4x4 grid in Supp. Fig. 6c, found {len(vc)}x{len(hr)}")

    bx, ax = np.polyfit(vc, x_labels, 1)
    by, ay = np.polyfit(hr, y_labels, 1)
    res_x = float(np.max(np.abs(np.array(vc) - (np.array(x_labels) - ax) / bx)))
    res_y = float(np.max(np.abs(np.array(hr) - (np.array(y_labels) - ay) / by)))
    if res_x > 4.0 or res_y > 4.0:
        raise RuntimeError(f"Grid calibration residual too high: rx={res_x:.2f}px, ry={res_y:.2f}px")

    r0, r1 = int(hr[0]) + 4, int(hr[-1]) - 4
    c0, c1 = int(vc[0]) + 4, int(vc[-1]) - 4
    sub = img[r0:r1, c0:c1]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]

    # 1) Blue filled circles (1-M, 1-H): medium-light blue interior + border, opened by 15x15
    blue_all = (B > 150.0) & (B - R > 40.0)
    circles_m = ndimage.binary_opening(blue_all, structure=np.ones((15, 15), dtype=bool))
    pts = []
    lab_c, nc = ndimage.label(circles_m)
    for k in range(1, nc + 1):
        rr, cc = np.nonzero(lab_c == k)
        if len(rr) < 700:
            continue
        c_abs = 0.5 * (float(cc.min()) + float(cc.max())) + 0.5 + c0
        r_abs = 0.5 * (float(rr.min()) + float(rr.max())) + 0.5 + r0
        pts.append((c_abs, r_abs, "circle", len(rr)))

    # 2) Dark-blue thick 'X' markers (a3-1, a4-1, a5-2, a6-1, a6-3): thick strokes (~10-12px) crossing at center (~18x18 px solid core!)
    # Notice italic blue text has thin strokes (~6px) and smaller color/saturation, whereas the 5 'X' markers are bold dark steel-blue (B - R > 50, G < 155, R < 110) with symmetric X-shaped ~55x55px connected components!
    dark_blue = (B - R > 45.0) & (B > 140.0) & (R < 115.0) & (G < 155.0)
    # Remove circle borders around already-found circle centers
    for c_abs, r_abs, _, _ in list(pts):
        rr_y, cc_x = np.ogrid[: dark_blue.shape[0], : dark_blue.shape[1]]
        dark_blue[(rr_y - (r_abs - r0)) ** 2 + (cc_x - (c_abs - c0)) ** 2 <= 42**2] = False
    # Each X marker has a solid crossing core of >= 11x11 px; opening by 9x9 isolates the 5 X crossing cores and kills all thin italic text!
    x_cores = ndimage.binary_opening(dark_blue, structure=np.ones((9, 9), dtype=bool))
    lab_x, nx = ndimage.label(x_cores)
    for k in range(1, nx + 1):
        rr, cc = np.nonzero(lab_x == k)
        if len(rr) < 120:
            continue
        c_abs = float(cc.mean()) + 0.5 + c0
        r_abs = float(rr.mean()) + 0.5 + r0
        pts.append((c_abs, r_abs, "X-cross", len(rr)))

    out = []
    for c_abs, r_abs, kind, area in pts:
        x_val = ax + bx * c_abs
        y_val = ay + by * r_abs
        err_px = 2.5 + res_y
        err_y = abs(by) * err_px
        out.append({
            "x": x_val,
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "px": [c_abs, r_abs],
            "px_lo": r_abs - err_px,
            "px_hi": r_abs + err_px,
            "note": f"{kind} marker (area={area}px, rx={res_x:.2f}px, ry={res_y:.2f}px)",
        })
    out.sort(key=lambda p: (round(p["x"], 1), p["y"]))
    if len(out) != 7:
        raise RuntimeError(f"Expected 7 markers in Supp. Fig. 6c, found {len(out)}: {[(round(p['x'],1), round(p['y'],1), p['note']) for p in out]}")
    return out
