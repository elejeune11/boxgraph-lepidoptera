"""Calibrated reader for Fig. 5(b) bubble radius distribution P(R) at t = 32768 (blue curve, p. 10, 1800 dpi)."""
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
    xt_inner = [c + ic0 + 8 for c, w in _clusters(np.where(top_band > 0.70)[0]) if w <= 14]
    x_val = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], dtype=float)
    x_px = np.array(
        [c0] + [min(xt_inner, key=lambda p: abs(p - (c0 + k * (c1 - c0) / 9.0))) for k in range(1, 9)] + [c1],
        dtype=float,
    )
    bx, ax = np.polyfit(x_px, x_val, 1)
    x_res = float(np.max(np.abs(x_px - (x_val - ax) / bx)))

    left_band = dark[ir0 + 8 : ir1 - 8, ic0 + 6 : ic0 + 18].mean(axis=1)
    yt_inner = [r + ir0 + 8 for r, w in _clusters(np.where(left_band > 0.70)[0]) if w <= 14]
    y_val = np.array([0.25, 0.20, 0.15, 0.10, 0.05, 0.0], dtype=float)
    y_px = np.array(
        [r0] + [min(yt_inner, key=lambda p: abs(p - (r0 + k * (r1 - r0) / 5.0))) for k in range(1, 5)] + [r1],
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

    def x2c(v): return int(round((v - ax) / bx))
    def y2r(v): return int(round((v - ay) / by))

    mask = np.zeros(R.shape, dtype=bool)
    mask[int(r0) + 6 : int(r1) - 3, int(c0) + 6 : int(c1) - 6] = True
    mask[: y2r(0.18), x2c(5.7) :] = False
    curve = mask & (B - R > 70) & (B - G > 70)

    r_lo, r_hi = int(r0) + 6, int(r1) - 2
    out = []
    for x in at:
        col_f = (float(x) - ax) / bx
        c = int(round(col_f))
        rs = np.where(curve[r_lo:r_hi, c])[0] + r_lo
        if len(rs) == 0:
            ys, xs = np.nonzero(curve[r_lo:r_hi, c - 8 : c + 9])
            if len(ys) == 0:
                out.append({"x": float(x), "y": None, "found": False, "note": "blue curve not found"})
                continue
            ys_abs = ys + r_lo
            xs_abs = xs + (c - 8)
            if len(np.unique(xs_abs)) >= 2:
                m, b = np.polyfit(xs_abs, ys_abs, 1)
                r_mid = float(m * col_f + b) + 0.5
                slope_px = abs(float(m))
            else:
                r_mid = float(ys_abs.mean()) + 0.5
                slope_px = 0.0
            half_thick = 6.0
        else:
            r_mid = float(rs.mean()) + 0.5
            half_thick = 0.5 * (float(rs.max() - rs.min()) + 1.0)
            rsl = np.where(curve[r_lo:r_hi, c - 6])[0] + r_lo
            rsr = np.where(curve[r_lo:r_hi, c + 6])[0] + r_lo
            slope_px = abs((rsr.mean() - rsl.mean()) / 12.0) if (len(rsl) > 0 and len(rsr) > 0) else 0.0
        err_px = half_thick + 1.0 + y_res + slope_px * (0.5 + x_res)
        y = max(0.0, float(ay + by * r_mid))
        v1 = max(0.0, float(ay + by * (r_mid - err_px)))
        v2 = max(0.0, float(ay + by * (r_mid + err_px)))
        out.append({
            "x": float(x),
            "y": y,
            "lo": min(v1, v2),
            "hi": max(v1, v2),
            "px": [float(col_f), r_mid],
            "px_lo": r_mid - err_px,
            "px_hi": r_mid + err_px,
            "note": f"blue t=32768 curve; cal residual ({x_res:.2f}, {y_res:.2f}) px",
        })
    return out
