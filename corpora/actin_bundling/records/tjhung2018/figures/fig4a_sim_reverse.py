"""Calibrated reader for Fig. 4(a) red square simulation markers (reverse Ostwald, p. 8, 1200 dpi)."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import find_objects, label


def _clusters(idxs: np.ndarray, gap: int = 3) -> list[tuple[float, int]]:
    if len(idxs) == 0:
        return []
    splits = np.where(np.diff(idxs) > gap)[0] + 1
    return [(float(g.mean()) + 0.5, int(len(g))) for g in np.split(idxs, splits)]


def _calibrate(img: np.ndarray):
    gray = img.astype(np.float32).mean(axis=2)
    dark = gray < 120.0
    r_lines = [c for c, _ in _clusters(np.where(dark.mean(axis=1) > 0.40)[0])]
    c_lines = [c for c, _ in _clusters(np.where(dark.mean(axis=0) > 0.40)[0])]
    r0, r1 = r_lines[0], r_lines[-1]
    c0, c1 = c_lines[0], c_lines[-1]
    ir0, ir1, ic0, ic1 = int(round(r0)), int(round(r1)), int(round(c0)), int(round(c1))

    top_band = dark[ir0 + 6 : ir0 + 18, ic0 + 8 : ic1 - 8].mean(axis=0)
    xt_inner = [c + ic0 + 8 for c, w in _clusters(np.where(top_band > 0.70)[0]) if w <= 12]
    x_px = np.array([c0] + xt_inner + [c1], dtype=float)
    x_val = np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0], dtype=float)
    bx, ax = np.polyfit(x_px, x_val, 1)
    x_res = float(np.max(np.abs(x_px - (x_val - ax) / bx)))

    left_band = dark[ir0 + 8 : ir1 - 8, ic0 + 6 : ic0 + 18].mean(axis=1)
    yt_cands = [r + ir0 + 8 for r, w in _clusters(np.where(left_band > 0.70)[0]) if w <= 12]
    y_val = np.array([2.0, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0], dtype=float)
    y_px_list = [r0]
    for k in range(1, 8):
        est = r0 + k * (r1 - r0) / 8.0
        best = min(yt_cands, key=lambda p: abs(p - est))
        y_px_list.append(best)
    y_px_list.append(r1)
    y_px = np.array(y_px_list, dtype=float)
    by, ay = np.polyfit(y_px, y_val, 1)
    y_res = float(np.max(np.abs(y_px - (y_val - ay) / by)))
    return (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1) = _calibrate(image)
    R = image[:, :, 0].astype(np.int16)
    G = image[:, :, 1].astype(np.int16)
    B = image[:, :, 2].astype(np.int16)

    def x2c(v): return int(round((v - ax) / bx))
    def y2r(v): return int(round((v - ay) / by))

    valid = np.zeros(R.shape, dtype=bool)
    # Include right frame edge so the half-clipped red square centered at (zeta=4, lam=1) is captured!
    valid[int(r0) + 5 : int(r1) - 5, int(c0) + 5 : int(c1) + 4] = True
    # Exclude insets and the large red region letters A, B, C
    valid[y2r(1.48) : y2r(0.30), x2c(1.70) : x2c(3.88)] = False
    valid[y2r(-0.56) : y2r(-1.80), x2c(0.55) : x2c(3.25)] = False
    valid[y2r(-0.80) : y2r(-2.00), x2c(-3.10) : x2c(-0.70)] = False

    red = valid & (R - G > 80) & (R - B > 80)
    lbl, _ = label(red)
    objs = find_objects(lbl)
    pts = []
    # Measure full square size from interior squares (~49x49 px) so border-clipped squares are centered accurately
    for k, sl in enumerate(objs, 1):
        if sl is None:
            continue
        sub = lbl[sl] == k
        area = int(sub.sum())
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        # Reject letters A, B, C (which are >100 px tall) and tiny noise
        if not (40 <= h <= 60 and 15 <= w <= 60 and area >= 600):
            continue
        ry, rx = np.nonzero(sub)
        rc = sl[0].start + float(ry.mean()) + 0.5
        if sl[1].stop >= int(round(c1)) - 5 and w < 0.7 * h:
            # Square is clipped by the right axis border at c1: center lies at sl[1].start + h/2 (on border c1)
            cc = float(sl[1].start) + 0.5 * h
            note = f"right-border-clipped red square ({w}x{h} px, area={area}); cal res ({x_res:.2f}, {y_res:.2f}) px"
        else:
            cc = sl[1].start + float(rx.mean()) + 0.5
            note = f"red square marker ({w}x{h} px, area={area}); cal res ({x_res:.2f}, {y_res:.2f}) px"
        x = float(ax + bx * cc)
        y = float(ay + by * rc)
        err_px = 3.0 + y_res
        v1, v2 = float(ay + by * (rc - err_px)), float(ay + by * (rc + err_px))
        pts.append({
            "x": x,
            "y": y,
            "lo": min(v1, v2),
            "hi": max(v1, v2),
            "px": [cc, rc],
            "px_lo": rc - err_px,
            "px_hi": rc + err_px,
            "note": note,
        })
    return sorted(pts, key=lambda d: d["x"])
