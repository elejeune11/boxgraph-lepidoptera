"""Calibrated reader for Fig. 6 Normal Ostwald (y=32) relative signed vapor current J_x(x)/|J_x(68)| (p. 16, 1500 dpi)."""
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

    bot_band = dark[ir1 - 18 : ir1 - 6, ic0 + 8 : ic1 - 8].mean(axis=0)
    xt_inner = [c + ic0 + 8 for c, w in _clusters(np.where(bot_band > 0.70)[0]) if w <= 12]
    x_px = np.array([c0] + xt_inner, dtype=float)
    x_val = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0], dtype=float)
    bx, ax = np.polyfit(x_px, x_val, 1)
    x_res = float(np.max(np.abs(x_px - (x_val - ax) / bx)))

    right_band = dark[ir0 + 8 : ir1 - 8, ic1 - 18 : ic1 - 6].mean(axis=1)
    yt_inner = [r + ir0 + 8 for r, w in _clusters(np.where(right_band > 0.70)[0]) if w <= 12]
    y_px = np.array(yt_inner + [r1], dtype=float)
    y_val = np.array([60.0, 50.0, 40.0, 30.0, 20.0, 10.0, 0.0], dtype=float)
    by, ay = np.polyfit(y_px, y_val, 1)
    y_res = float(np.max(np.abs(y_px - (y_val - ay) / by)))
    return (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1) = _calibrate(image)
    R = image[:, :, 0].astype(np.int16)
    G = image[:, :, 1].astype(np.int16)
    B = image[:, :, 2].astype(np.int16)

    def y2r(v): return int(round((v - ay) / by))

    # Extract isolated horizontal vapor arrows along the midline Normal Ostwald row y ≈ 32.0
    r_top, r_bot = y2r(33.5), y2r(30.5)
    c_l, c_r = int(round(c0)) + 8, int(round(c1)) - 8
    red = (R[r_top:r_bot, c_l:c_r] - B[r_top:r_bot, c_l:c_r] > 60) & (
        R[r_top:r_bot, c_l:c_r] - G[r_top:r_bot, c_l:c_r] > 60
    )
    lbl, _ = label(red)
    objs = find_objects(lbl)
    arrows = []
    for k, sl in enumerate(objs, 1):
        if sl is None:
            continue
        sub = lbl[sl] == k
        area = int(sub.sum())
        w_px = sl[1].stop - sl[1].start
        h_px = sl[0].stop - sl[0].start
        if area < 400 or w_px < 40 or w_px > 180 or h_px > 42:
            continue
        ry, rx = np.nonzero(sub)
        cc = c_l + sl[1].start + float(rx.mean()) + 0.5
        rc = r_top + sl[0].start + float(ry.mean()) + 0.5
        x = float(ax + bx * cc)
        # Keep clean vapor-region arrows (inter-droplet gap x in [58, 86] and outer vapor x in [110, 126])
        if not ((58.0 <= x <= 86.0) or (110.0 <= x <= 126.0)):
            continue
        # Direction from arrowhead vertical spread (left 35% vs right 35% of arrow)
        ml = rx < 0.35 * w_px
        mr = rx > 0.65 * w_px
        h_l = int(np.ptp(ry[ml])) if ml.any() else 0
        h_r = int(np.ptp(ry[mr])) if mr.any() else 0
        sign = -1.0 if h_l > h_r else 1.0
        # Approximate arrow base x (tail) vs tip: shaft length ~ w_px - 6 px (stroke width)
        length = float(w_px) - 5.0
        arrows.append((x, sign, length, cc, rc, w_px))

    arrows.sort(key=lambda t: t[0])
    # Reference arrow closest to mid-gap x = 68.0
    ref = min(arrows, key=lambda t: abs(t[0] - 68.0))
    ref_len = ref[2]

    out = []
    for x, sign, length, cc, rc, w_px in arrows:
        val = sign * (length / ref_len)
        # Honest relative uncertainty: +/- 4 px length uncertainty + reference uncertainty
        rel_err = abs(val) * ((4.0 / length) ** 2 + (4.0 / ref_len) ** 2) ** 0.5
        out.append({
            "x": x,
            "y": val,
            "lo": val - rel_err,
            "hi": val + rel_err,
            "px": [cc, rc],
            "note": f"Normal Ostwald vapor arrow y≈32, w={w_px} px, dir={'left' if sign < 0 else 'right'}",
        })
    return out
