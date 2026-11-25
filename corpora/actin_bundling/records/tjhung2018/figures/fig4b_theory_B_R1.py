"""Calibrated reader for Fig. 4(b) Region B blue dashed theory curve R_1(t) (zeta=-4, lam=-1, p. 8, 1200 dpi)."""
from __future__ import annotations

import math
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

    top_major = dark[ir0 + 14 : ir0 + 26, ic0 + 8 : ic1 - 8].mean(axis=0)
    xt_major = [c + ic0 + 8 for c, w in _clusters(np.where(top_major > 0.70)[0]) if w <= 12]
    dec_est = (c1 - c0) / (3.0 + math.log10(2.0))
    x_px = np.array(
        [c0] + [min(xt_major, key=lambda p: abs(p - (c0 + k * dec_est))) for k in (1, 2, 3)],
        dtype=float,
    )
    x_u = np.array([2.0, 3.0, 4.0, 5.0], dtype=float)
    bx, ax = np.polyfit(x_px, x_u, 1)
    x_res = float(np.max(np.abs(x_px - (x_u - ax) / bx)))

    right_band = dark[ir0 + 8 : ir1 - 8, ic1 - 18 : ic1 - 6].mean(axis=1)
    yt_cands = [r + ir0 + 8 for r, w in _clusters(np.where(right_band > 0.70)[0]) if w <= 12]
    y_val = np.array([25.0, 20.0, 15.0, 10.0, 5.0, 0.0], dtype=float)
    y_px = np.array(
        [r0] + [min(yt_cands, key=lambda p: abs(p - (r0 + k * (r1 - r0) / 5.0))) for k in range(1, 5)] + [r1],
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

    def y2r(v): return int(round((v - ay) / by))

    mask = (B - R > 70) & (B - G > 70)
    # Upper Region B branch lies in R in [13.20, 18.0]
    r_top, r_bot = y2r(18.0), y2r(13.20)

    out = []
    for x in at:
        col_f = (math.log10(float(x)) - ax) / bx
        c = int(round(col_f))
        hw = 22
        ca, cb = max(int(c0) + 5, c - hw), min(int(c1) - 5, c + hw + 1)
        strip = mask[r_top:r_bot, ca:cb]
        proj = strip.sum(axis=1)
        rs = np.where(proj > 0)[0]
        if len(rs) == 0:
            out.append({"x": float(x), "y": None, "found": False, "note": "blue dashed R1(t) not found"})
            continue
        splits = np.where(np.diff(rs) > 5)[0] + 1
        bands = np.split(rs, splits)
        # In the upper branch [13.20, 18.0], solid R1 is topmost (bands[0]), dashed R1 is below it (bands[-1])!
        # At t > 110000 (after solid ends), only the converged dashed asymptote remains (bands[0] == bands[-1]).
        g = bands[-1]
        ys, xs = np.nonzero(strip[g[0] : g[-1] + 1, :])
        ys_abs = ys + g[0] + r_top
        xs_abs = xs + ca
        if len(np.unique(xs_abs)) >= 2:
            m, b = np.polyfit(xs_abs.astype(float), ys_abs.astype(float), 1)
            r_mid = float(m * col_f + b) + 0.5
            resid = np.abs(ys_abs - (m * xs_abs + b))
            half_thick = max(3.5, float(np.percentile(resid, 95)))
            slope_px = abs(float(m))
        else:
            r_mid = float(ys_abs.mean()) + 0.5
            half_thick = 4.0
            slope_px = 0.0
        err_px = half_thick + 1.0 + y_res + slope_px * (1.0 + x_res)
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
            "note": f"blue dashed R1(t); cal residual ({x_res:.2f}, {y_res:.2f}) px",
        })
    return out
