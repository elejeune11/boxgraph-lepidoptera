"""Calibrated reader for Fig. 2C (PDF p. 3 / printed p. 26572):
Binned median and 10th-90th percentile vertical scatter of lipid domain excess perimeter
L_ex = L/L_0 - 1 vs. normalized actin density rho/rho_0 (n = 270 open black square samples).
Includes sampling directly at the four vertical pink dashed lines (panels 1-4 of Fig. 2B),
located from the image at rho/rho_0 = 0.048, 0.372, 0.674, and 0.976.
"""

from __future__ import annotations

import numpy as np

from boxgraph import digitize as D


def _group(indices: np.ndarray, gap: int = 2) -> list[list[int]]:
    out: list[list[int]] = []
    for x in indices:
        xi = int(x)
        if out and xi <= out[-1][-1] + gap:
            out[-1].append(xi)
        else:
            out.append([xi])
    return out


def _calibrate(image: np.ndarray):
    H, W = image.shape[:2]
    r, g, b = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    black = (mx < 100) & ((mx - mn) < 25)

    # Bottom x-axis spine and left y-axis spine
    bot_g = _group(np.where(black.sum(axis=1) > 0.5 * W)[0], gap=4)[-1]
    left_g = _group(np.where(black.sum(axis=0) > 0.5 * H)[0], gap=4)[0]
    r_x0, r_x1 = bot_g[0], bot_g[-1]
    c_y0, c_y1 = left_g[0], left_g[-1]
    r_zero = (r_x0 + r_x1) / 2.0 + 0.5
    c_zero = (c_y0 + c_y1) / 2.0 + 0.5

    # X-axis ticks pointing upward from bottom spine: 0.2, 0.4, 0.6, 0.8, 1.0
    xt_raw = []
    for c in range(c_y1 + 10, W - 5):
        n = 0
        while r_x0 - 1 - n >= 0 and black[r_x0 - 1 - n, c]:
            n += 1
        if 12 <= n <= 40:
            xt_raw.append(c)
    xt_groups = [grp for grp in _group(np.array(xt_raw), gap=1) if 4 <= len(grp) <= 18]
    xt_pos = sorted((grp[0] + grp[-1]) / 2.0 + 0.5 for grp in xt_groups)
    x_ticks = D.regular_ticks([c_zero] + xt_pos, tol=3.0)
    if len(x_ticks) != 6:
        raise D.CalibrationError(f"Fig. 2C: expected 6 x-ticks (0..1 step 0.2), got {x_ticks}")
    xa = D.Axis("linear", list(zip(x_ticks, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])))

    # Y-axis ticks pointing rightward from left spine: 1.0, 0.5 (length ~21 px, height ~9 px;
    # exclude data squares touching the spine which are ~48 px long and ~18 px thick)
    yt_raw = []
    for r_idx in range(100, r_x0 - 20):
        n = 0
        while c_y1 + 1 + n < W and black[r_idx, c_y1 + 1 + n]:
            n += 1
        if 14 <= n <= 32:
            yt_raw.append(r_idx)
    yt_groups = [grp for grp in _group(np.array(yt_raw), gap=1) if 4 <= len(grp) <= 14]
    yt_pos = sorted((grp[0] + grp[-1]) / 2.0 + 0.5 for grp in yt_groups)
    y_ticks = D.regular_ticks(yt_pos + [r_zero], tol=3.0)
    if len(y_ticks) != 3:
        raise D.CalibrationError(f"Fig. 2C: expected 3 y-ticks (1.0, 0.5, 0.0), got {y_ticks}")
    ya = D.Axis("linear", list(zip(y_ticks, [1.0, 0.5, 0.0])))

    return xa, ya, black, (c_y1, r_x0), yt_pos


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, black, (c_y1, r_x0), yt_pos = _calibrate(image)
    H, W = image.shape[:2]

    # Restrict to plot interior strictly inside spines, clearing the two interior y-axis ticks
    mask = np.zeros((H, W), dtype=bool)
    c_right_lim = min(W, int(round(xa.to_pixel(1.02))))
    mask[80 : r_x0 - 25, c_y1 + 2 : c_right_lim] = black[80 : r_x0 - 25, c_y1 + 2 : c_right_lim]
    for ry in yt_pos:
        ri = int(round(ry))
        mask[max(0, ri - 12) : min(H, ri + 13), c_y1 : c_y1 + 30] = False

    # Bin half-width +-0.03 in normalized actin density rho/rho_0 (~65 px)
    hw = int(round(0.03 * abs(1.0 / xa.b)))
    base_err_px = 2.0 + max(ya.residual_px, 0.5)

    out = []
    for xv in at:
        pc = float(xa.to_pixel(xv))
        c0 = int(round(pc))
        c_lo = max(c_y1 + 2, c0 - hw)
        c_hi = min(c_right_lim, c0 + hw + 1)
        rr = np.where(mask[:, c_lo:c_hi])[0]
        if len(rr) < 50:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no scatter squares near rho/rho_0={xv:g}"})
            continue
        r_med = float(np.median(rr)) + 0.5
        r_p10 = float(np.percentile(rr, 10)) + 0.5
        r_p90 = float(np.percentile(rr, 90)) + 0.5
        r_hi_edge = min(r_p10, r_med) - base_err_px
        r_lo_edge = max(r_p90, r_med) + base_err_px
        y_val = float(ya.to_value(r_med))
        a, b = ya.to_value(r_hi_edge), ya.to_value(r_lo_edge)
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [pc, r_med],
            "px_lo": float(r_hi_edge),
            "px_hi": float(r_lo_edge),
            "note": (
                f"binned median of scatter (bin +-0.03; 10-90 pct spread); "
                f"cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px"
            ),
        })
    return out
