"""Calibrated reader for Fig. 6 Reverse Ostwald (bottom half, y < 30) blue droplet interface contours y(x) (p. 16, 1500 dpi)."""
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
    xt_inner = [c + ic0 + 8 for c, w in _clusters(np.where(bot_band > 0.70)[0]) if w <= 12]
    assert len(xt_inner) == 6, f"expected 6 interior x-ticks (20..120), got {len(xt_inner)}"
    x_px = np.array([c0] + xt_inner, dtype=float)
    x_val = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0], dtype=float)
    bx, ax = np.polyfit(x_px, x_val, 1)
    x_res = float(np.max(np.abs(x_px - (x_val - ax) / bx)))

    right_band = dark[ir0 + 8 : ir1 - 8, ic1 - 18 : ic1 - 6].mean(axis=1)
    yt_inner = [r + ir0 + 8 for r, w in _clusters(np.where(right_band > 0.70)[0]) if w <= 12]
    assert len(yt_inner) == 6, f"expected 6 interior right y-ticks (60..10), got {len(yt_inner)}"
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

    blue = (B - R > 60) & (B - G > 40)
    r_midline = y2r(29.7)
    r_bot = int(round(r1)) - 5

    out = []
    for x in at:
        col_f = (float(x) - ax) / bx
        c = int(round(col_f))
        rs = np.where(blue[r_midline:r_bot, c])[0] + r_midline
        if len(rs) == 0:
            rs = np.where(blue[r_midline:r_bot, c - 4 : c + 5].any(axis=1))[0] + r_midline
        if len(rs) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": "bottom blue contour not found"})
            continue
        r_mid = float(rs.mean()) + 0.5
        half_thick = 0.5 * (float(rs.max() - rs.min()) + 1.0)
        slope_px = 0.0
        rsl = np.where(blue[r_midline:r_bot, c - 5])[0] + r_midline
        rsr = np.where(blue[r_midline:r_bot, c + 5])[0] + r_midline
        if len(rsl) > 0 and len(rsr) > 0:
            slope_px = abs((rsr.mean() - rsl.mean()) / 10.0)
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
            "note": f"Reverse Ostwald bottom droplet contour; thickness {len(rs)} px; cal res ({x_res:.2f}, {y_res:.2f}) px",
        })
    return out
