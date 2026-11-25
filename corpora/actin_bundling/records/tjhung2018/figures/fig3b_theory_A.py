"""Calibrated reader for Fig. 3(b) Region A analytical curve phi_+(R) (red, zeta=-1, lam=0.5, p. 6, 1200 dpi)."""
from __future__ import annotations

import numpy as np


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
    xt_cands = [c + ic0 + 8 for c, w in _clusters(np.where(bot_band > 0.70)[0]) if w <= 12]
    x_val = np.array([0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0], dtype=float)
    x_px = np.array(
        [c0] + [min(xt_cands, key=lambda p: abs(p - (c0 + k * (c1 - c0) / 8.0))) for k in range(1, 8)] + [c1],
        dtype=float,
    )
    bx, ax = np.polyfit(x_px, x_val, 1)
    x_res = float(np.max(np.abs(x_px - (x_val - ax) / bx)))

    left_band = dark[ir0 + 8 : ir1 - 8, ic0 + 6 : ic0 + 18].mean(axis=1)
    yt_cands = [r + ir0 + 8 for r, w in _clusters(np.where(left_band > 0.70)[0]) if w <= 12]
    y_val = np.array([1.3, 1.2, 1.1, 1.0], dtype=float)
    y_px = np.array(
        [r0] + [min(yt_cands, key=lambda p: abs(p - (r0 + k * (r1 - r0) / 3.0))) for k in range(1, 3)] + [r1],
        dtype=float,
    )
    by, ay = np.polyfit(y_px, y_val, 1)
    y_res = float(np.max(np.abs(y_px - (y_val - ay) / by)))
    return (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    (ax, bx, x_res), (ay, by, y_res), (c0, r0, c1, r1) = _calibrate(image)
    R = image[:, :, 0].astype(np.int16)
    G = image[:, :, 1].astype(np.int16)
    B = image[:, :, 2].astype(np.int16)
    mask = (R - G > 70) & (R - B > 70)

    r_lo, r_hi = int(round(r0)) + 5, int(round(r1)) - 5
    out = []
    for x in at:
        col_f = (float(x) - ax) / bx
        c = int(round(col_f))
        rs = np.where(mask[r_lo:r_hi, c])[0] + r_lo
        if len(rs) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": "red curve not found"})
            continue
        splits = np.where(np.diff(rs) > 6)[0] + 1
        # Red theory curve lies at/above any red '+' markers
        g = np.split(rs, splits)[0]
        r_mid = float(g.mean()) + 0.5
        half_thick = 0.5 * (float(g.max() - g.min()) + 1.0)
        slope_px = 0.0
        c_l, c_r = max(int(c0) + 10, c - 6), min(int(c1) - 10, c + 6)
        rsl = np.where(mask[r_lo:r_hi, c_l])[0] + r_lo
        rsr = np.where(mask[r_lo:r_hi, c_r])[0] + r_lo
        if len(rsl) > 0 and len(rsr) > 0 and c_r > c_l:
            gl = np.split(rsl, np.where(np.diff(rsl) > 6)[0] + 1)[0]
            gr = np.split(rsr, np.where(np.diff(rsr) > 6)[0] + 1)[0]
            slope_px = abs((gr.mean() - gl.mean()) / (c_r - c_l))
        err_px = half_thick + 0.5 + y_res + slope_px * (0.5 + x_res)
        y = float(ay + by * r_mid)
        v1, v2 = float(ay + by * (r_mid - err_px)), float(ay + by * (r_mid + err_px))
        out.append({
            "x": float(x),
            "y": y,
            "lo": min(v1, v2),
            "hi": max(v1, v2),
            "px": [float(col_f), r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"red curve; thickness {len(g)} px; cal residual ({x_res:.2f}, {y_res:.2f}) px",
        })
    return out
