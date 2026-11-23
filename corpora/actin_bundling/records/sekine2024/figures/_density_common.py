"""Shared calibration and kernel-density curve reader for Stage 15 (red) & Stage 16 (teal) plots (Fig. 1f, Supp. Fig. 1b,c)."""

from __future__ import annotations

import numpy as np


def _split_runs(idx: np.ndarray, gap: int = 5) -> list[np.ndarray]:
    if len(idx) == 0:
        return []
    cuts = np.where(np.diff(idx) > gap)[0] + 1
    return [g for g in np.split(idx, cuts) if len(g) > 0]


def read_density_curve(
    image: np.ndarray,
    at: list[float],
    x_labels: list[float],
    y_labels_top_to_bottom: list[float],
    curve: str = "st15_red",
    tick_band_x: tuple[int, int] = (10, 35),
    tick_band_y: tuple[int, int] = (10, 35),
):
    """Calibrate detached L-axes and read either the vivid red (St 15) or vivid teal (St 16) kernel density curve."""
    img = image.astype(float)
    H, W, _ = img.shape
    gray = img.mean(axis=2)
    dark = gray < 150.0

    hr = _split_runs(np.where(dark.sum(axis=1) > W * 0.40)[0], gap=10)[-1]
    vc = _split_runs(np.where(dark.sum(axis=0) > H * 0.40)[0], gap=10)[0]

    xp = (255.0 - gray[int(hr.max()) + tick_band_x[0] : min(H, int(hr.max()) + tick_band_x[1]), :]).mean(axis=0)
    xp[: int(vc.max()) + 2] = 0.0
    xt = [float(((g + 0.5) * xp[g]).sum() / xp[g].sum()) for g in _split_runs(np.where(xp > 45.0)[0], gap=5)]
    if len(xt) != len(x_labels):
        raise RuntimeError(f"Expected {len(x_labels)} x-ticks, found {len(xt)}: {xt}")

    c1_yt = max(0, int(vc.min()) - tick_band_y[0])
    c0_yt = max(0, int(vc.min()) - tick_band_y[1])
    yp = (255.0 - gray[:, c0_yt:c1_yt]).mean(axis=1)
    yp[int(hr.min()) + 10 :] = 0.0
    yt = [float(((g + 0.5) * yp[g]).sum() / yp[g].sum()) for g in _split_runs(np.where(yp > 45.0)[0], gap=5)]
    if len(yt) != len(y_labels_top_to_bottom):
        raise RuntimeError(f"Expected {len(y_labels_top_to_bottom)} y-ticks, found {len(yt)}: {yt}")

    bx, ax = np.polyfit(xt, x_labels, 1)
    by, ay = np.polyfit(yt, y_labels_top_to_bottom, 1)
    res_x = float(np.max(np.abs(np.array(xt) - (np.array(x_labels) - ax) / bx)))
    res_y = float(np.max(np.abs(np.array(yt) - (np.array(y_labels_top_to_bottom) - ay) / by)))
    if res_x > 2.5 or res_y > 2.5:
        raise RuntimeError(f"Calibration residual too high: rx={res_x:.2f}px, ry={res_y:.2f}px")

    r_top = max(10, int(yt[0]) - 60)
    r_bot = min(int(hr.min()) - 4, int(yt[-1]) + 20)
    c_left = int(vc.max()) + 4
    c_right = W - 10
    sub = img[r_top:r_bot, c_left:c_right]
    R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]

    if curve == "st15_red":
        # Vivid orange-red density curve (distinct from pale pink histogram fill/border)
        mask = (R - G > 90.0) & (R - B > 100.0) & (G < 130.0)
    elif curve == "st16_teal":
        # Vivid teal density curve (distinct from pale cyan histogram fill/border)
        mask = (G - R > 60.0) & (B - R > 60.0) & (R < 90.0)
    else:
        raise ValueError(f"Unknown curve: {curve}")

    # Exclude legend box region on the upper right (above y = 0.45*y_max and right of 55% of x_max)
    c_leg = int(round((0.52 * max(x_labels) - ax) / bx)) - c_left
    r_leg = int(round((0.35 * max(y_labels_top_to_bottom) - ay) / by)) - r_top
    if c_leg > 0 and r_leg > 0:
        mask[:r_leg, c_leg:] = False

    W_s = mask.shape[1]
    out = []
    for xv in at:
        xv = float(xv)
        col_f = (xv - ax) / bx
        col_s = int(round(col_f)) - c_left
        if not (0 <= col_s < W_s):
            out.append({"x": xv, "y": None, "note": "column out of bounds"})
            continue
        rr = np.nonzero(mask[:, col_s])[0]
        if len(rr) == 0:
            # Try ±2 px window in case of slight antialiasing break
            rr, _ = np.nonzero(mask[:, max(0, col_s - 2) : min(W_s, col_s + 3)])
        if len(rr) == 0:
            out.append({"x": xv, "y": None, "note": f"{curve} curve not present at x={xv}"})
            continue
        runs = _split_runs(np.unique(rr), gap=8)
        best = max(runs, key=len)
        rc_abs = 0.5 * (float(best[0]) + float(best[-1])) + 0.5 + r_top
        half_w = 0.5 * float(best[-1] - best[0] + 1)
        y_val = ay + by * rc_abs
        err_px = max(2.0, 0.5 * half_w) + res_y
        err_y = abs(by) * err_px
        out.append({
            "x": xv,
            "y": y_val,
            "lo": y_val - err_y,
            "hi": y_val + err_y,
            "px": [col_f, rc_abs],
            "px_lo": rc_abs - err_px,
            "px_hi": rc_abs + err_px,
            "note": f"{curve} density (span={int(best[-1]-best[0]+1)}px, rx={res_x:.2f}px, ry={res_y:.2f}px)",
        })
    return out
