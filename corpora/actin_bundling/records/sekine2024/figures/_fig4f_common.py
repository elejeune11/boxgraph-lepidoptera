"""Shared dual-axis reader for Fig. 4f (4-M orientation and coherency time courses, t = 0..40 s, PDF p. 7, 1200 dpi)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read_fig4f(image: np.ndarray, series: str):
    """Read all 41 time points (t = 0, 1, ..., 40 s) of either 'orientation' (green, 0-90 deg) or 'coherency' (blue, 0-0.3)."""
    img = image.astype(float)
    H, W, _ = img.shape
    g = img.mean(axis=2)

    is_grid = (g < 235.0) & (g > 150.0) & (np.abs(img[:, :, 0] - img[:, :, 2]) < 10.0)
    hr = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_grid.sum(axis=1) > W * 0.50)[0], gap=10)]
    vc = [float(gr.mean()) + 0.5 for gr in _split_runs(np.where(is_grid.sum(axis=0) > H * 0.40)[0], gap=10)]
    if len(hr) != 4 or len(vc) != 2:
        raise RuntimeError(f"Expected 4 horizontal gridlines and 2 vertical borders, found {len(hr)}x{len(vc)}")

    if series == "orientation":
        y_labels = [90.0, 60.0, 30.0, 0.0]
    elif series == "coherency":
        y_labels = [0.3, 0.2, 0.1, 0.0]
    else:
        raise ValueError(f"Unknown series: {series}")

    by, ay = np.polyfit(hr, y_labels, 1)
    res_y = float(np.max(np.abs(np.array(hr) - (np.array(y_labels) - ay) / by)))
    if res_y > 2.5:
        raise RuntimeError(f"Y calibration residual too high: ry={res_y:.2f}px")

    # Excel category axis: 41 equally spaced bins across [vc[0], vc[1]] with tick k*10 at left edge of bin k*10
    c_left, c_right = vc[0], vc[1]
    dx = (c_right - c_left) / 41.0

    r0 = max(10, int(hr[0]) - 28)
    r1 = min(H - 10, int(hr[-1]) + 28)
    sub = img[r0:r1, :]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]

    if series == "orientation":
        mask = (G - R > 15.0) & (G - B > 25.0)
    else:
        mask = (B - R > 25.0) & (B > 140.0)

    # First detect clean isolated circle centers via 17x17 opening (reject tall blobs > 42px where sharp V-legs merged with the circle)
    opened = ndimage.binary_opening(mask, structure=np.ones((17, 17), dtype=bool))
    lab, n = ndimage.label(opened)
    detected = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if len(rr) > 150 and (rr.max() - rr.min() + 1) <= 42 and (cc.max() - cc.min() + 1) <= 45:
            detected.append((float(cc.mean()) + 0.5, float(rr.mean()) + 0.5 + r0))

    # Match each of the 41 integer seconds t = 0..40 to its circle center (or column fallback if overlapping)
    out = []
    for t in range(41):
        c_exp = c_left + (t + 0.5) * dx
        near = [p for p in detected if abs(p[0] - c_exp) <= 14.0]
        if near:
            p = min(near, key=lambda z: abs(z[0] - c_exp))
            c_use, r_use = p[0], p[1]
            err_px = 2.0 + res_y
            src = "circle"
        else:
            # Fallback for overlapping or steeply joined vertex (e.g. t=1 where green and blue meet near 0)
            ci = int(round(c_exp))
            rr = np.nonzero(mask[:, max(0, ci - 4) : min(W, ci + 5)].any(axis=1))[0]
            if len(rr) == 0:
                out.append({"x": float(t), "y": None, "note": f"no {series} pixels at t={t}"})
                continue
            runs = _split_runs(rr, gap=6)
            best = max(runs, key=len)
            # If steep V-trough or peak, pick trough/peak cap using circle radius ~17px
            if len(best) > 45:
                # Check neighbors to tell peak vs trough
                c_p = int(round(c_exp - dx)); c_n = int(round(c_exp + dx))
                rr_p = np.nonzero(mask[:, max(0, c_p - 3) : min(W, c_p + 4)].any(axis=1))[0]
                rr_n = np.nonzero(mask[:, max(0, c_n - 3) : min(W, c_n + 4)].any(axis=1))[0]
                mp = float(rr_p.mean()) if len(rr_p) else float(best.mean())
                mn = float(rr_n.mean()) if len(rr_n) else float(best.mean())
                if mp > best[0] + 18 and mn > best[0] + 18:
                    r_use = float(best[0]) + 17.0 + r0
                elif mp < best[-1] - 18 and mn < best[-1] - 18:
                    r_use = float(best[-1]) - 17.0 + r0
                else:
                    r_use = 0.5 * (float(best[0]) + float(best[-1])) + 0.5 + r0
            else:
                r_use = 0.5 * (float(best[0]) + float(best[-1])) + 0.5 + r0
            c_use = c_exp
            err_px = 4.0 + res_y
            src = "col-fallback"

        y_val = ay + by * r_use
        err_y = abs(by) * err_px
        out.append({
            "x": float(t),
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "px": [c_use, r_use],
            "px_lo": r_use - err_px,
            "px_hi": r_use + err_px,
            "note": f"4-M {series} ({src}, ry={res_y:.2f}px)",
        })
    return out
