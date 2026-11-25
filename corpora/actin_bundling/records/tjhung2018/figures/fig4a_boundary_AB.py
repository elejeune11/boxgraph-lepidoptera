"""Calibrated reader for Fig. 4(a) Region A/B analytical boundary lambda*(zeta) (p. 8, 1200 dpi)."""
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

    top_band = dark[ir0 + 6 : ir0 + 18, ic0 + 8 : ic1 - 8].mean(axis=0)
    xt_inner = [c + ic0 + 8 for c, w in _clusters(np.where(top_band > 0.70)[0]) if w <= 12]
    x_px = np.array([c0] + xt_inner + [c1], dtype=float)
    x_val = np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0], dtype=float)
    assert len(x_px) == len(x_val), f"expected 9 x-ticks, got {len(x_px)}"
    bx, ax = np.polyfit(x_px, x_val, 1)
    x_res = float(np.max(np.abs(x_px - (x_val - ax) / bx)))

    # Match 7 inner y-ticks against expected positions from [r0, r1]
    left_band = dark[ir0 + 8 : ir1 - 8, ic0 + 6 : ic0 + 18].mean(axis=1)
    yt_cands = [r + ir0 + 8 for r, w in _clusters(np.where(left_band > 0.70)[0]) if w <= 12]
    y_val = np.array([2.0, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0], dtype=float)
    y_px_list = [r0]
    for k in range(1, 8):
        est = r0 + k * (r1 - r0) / 8.0
        best = min(yt_cands, key=lambda p: abs(p - est))
        assert abs(best - est) < 8.0, f"missing y-tick near {est:.1f}"
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
    blk = (R < 90) & (G < 90) & (B < 90)

    r_zero = int(round((0.0 - ay) / by))
    r_bot = int(round(r1)) - 15
    # The Region A/B curve (zeta < 0, lam < 0) is the topmost black curve below lam = 0
    # (clear of the left tick at lam = 0 by starting 8 px below r_zero).
    r_start = r_zero + 8

    out = []
    for x in at:
        col_f = (float(x) - ax) / bx
        c = int(round(col_f))
        rs = np.where(blk[r_start:r_bot, c])[0] + r_start
        if len(rs) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": "black curve not found"})
            continue
        splits = np.where(np.diff(rs) > 5)[0] + 1
        g = np.split(rs, splits)[0]
        r_mid = float(g.mean()) + 0.5
        half_thick = 0.5 * (float(g.max() - g.min()) + 1.0)
        # Estimate slope across +/- 6 px
        slope_px = 0.0
        c_l, c_r = max(int(c0) + 10, c - 6), min(int(c1) - 10, c + 6)
        rsl = np.where(blk[r_start:r_bot, c_l])[0] + r_start
        rsr = np.where(blk[r_start:r_bot, c_r])[0] + r_start
        if len(rsl) > 0 and len(rsr) > 0 and c_r > c_l:
            gl = np.split(rsl, np.where(np.diff(rsl) > 5)[0] + 1)[0]
            gr = np.split(rsr, np.where(np.diff(rsr) > 5)[0] + 1)[0]
            slope_px = abs((gr.mean() - gl.mean()) / (c_r - c_l))
        err_px = half_thick + 0.5 + y_res + slope_px * (0.5 + x_res)
        y = float(ay + by * r_mid)
        v1 = float(ay + by * (r_mid - err_px))
        v2 = float(ay + by * (r_mid + err_px))
        out.append({
            "x": float(x),
            "y": y,
            "lo": min(v1, v2),
            "hi": max(v1, v2),
            "px": [float(col_f), r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"black A/B curve; thickness {len(g)} px; cal residual ({x_res:.2f}, {y_res:.2f}) px",
        })
    return out
