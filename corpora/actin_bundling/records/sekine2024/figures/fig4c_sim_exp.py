"""Digitize all 9 unoccluded markers (7 CGMD simulations 1-M, 2-M, 2-H, 3-M, 3-H, 4-M, 4-H + experimental St.15, St.16; 1-H is occluded behind 3-H) in Fig. 4c (PDF p. 7, 900 dpi)."""

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

    is_ng = (g < 238.0) & (g > 150.0) & (np.abs(img[:, :, 0] - img[:, :, 2]) < 12.0)
    vc = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_ng.sum(axis=0) > H * 0.50)[0], gap=10)]

    c_l, c_r = int(vc[0]) + 5, int(vc[-1]) - 5
    clean_cols = np.zeros(W, dtype=bool)
    clean_cols[c_l : c_l + 200] = True
    clean_cols[c_r - 200 : c_r] = True
    hr = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_ng[:, clean_cols].sum(axis=1) > 65)[0], gap=12)]

    x_labels = [30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    y_labels = [500.0, 450.0, 400.0, 350.0, 300.0, 250.0]
    if len(vc) != len(x_labels) or len(hr) != len(y_labels):
        raise RuntimeError(f"Expected 6x6 grid, found {len(vc)}x{len(hr)}")

    bx, ax = np.polyfit(vc, x_labels, 1)
    by, ay = np.polyfit(hr, y_labels, 1)
    res_x = float(np.max(np.abs(np.array(vc) - (np.array(x_labels) - ax) / bx)))
    res_y = float(np.max(np.abs(np.array(hr) - (np.array(y_labels) - ay) / by)))
    if res_x > 2.5 or res_y > 3.5:
        raise RuntimeError(f"Grid calibration residual too high: rx={res_x:.2f}px, ry={res_y:.2f}px")

    # Exclude top legend band (y > 460)
    r0 = int(round((460.0 - ay) / by))
    r1 = int(hr[-1]) - 2
    c0, c1 = int(vc[0]) + 2, int(vc[-1]) - 2
    sub = img[r0:r1, c0:c1]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    gsub = sub.mean(axis=2)

    # Disjoint marker-fill color masks (no fill_holes so closed digits '4','6' in text labels are never filled!)
    specs = [
        # Red diamonds (4-H, 4-M): pinkish-red fill + bright red border (B > 90 in salmon fill or R > 210, G < 100 in border)
        ("red_4", ((R > 205) & (R - G > 55) & (B > 80) & (G < 185)) | ((R > 210) & (G < 90) & (B < 90)), 13, 600),
        # Orange triangles (2-M, 2-H): orange border + light peach/orange interior (R > 210, 95 < G < 210, B < 140)
        ("orange_2", (R > 210) & (G > 95) & (G < 205) & (B < 145) & (R - B > 95) & (R - G > 25), 9, 450),
        # Yellow squares (3-H, 3-M): yellow interior + golden border
        ("yellow_3", (R > 180) & (G > 165) & (B < 170) & (G - B > 40) & (R - G < 55), 11, 400),
        # Blue circle (1-M; 1-H is occluded directly behind 3-H)
        ("blue_1", (B > 175) & (B - R > 35) & (G > 125), 13, 700),
        # Gray experimental markers (St.15 circle, St.16 diamond): neutral gray fill (excluding thin dark text strokes)
        ("gray_exp", (gsub < 205) & (gsub > 105) & (np.abs(R - B) < 12) & (np.abs(R - G) < 12), 13, 600),
    ]

    pts = []
    for label_kind, mask, ksz, min_area in specs:
        # Close tiny 3px gaps inside markers only (never fill large letter loops)
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
            # Deduplicate if an adjacent color mask overlapped the same marker center within 15 px
            if any(abs(p["px"][0] - c_abs) < 15.0 and abs(p["px"][1] - r_abs) < 15.0 for p in pts):
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
    if len(pts) != 9:
        raise RuntimeError(f"Expected 9 unoccluded markers in Fig. 4c, found {len(pts)}: {[(round(p['x'],1), round(p['y'],1), p['note']) for p in pts]}")
    return pts
