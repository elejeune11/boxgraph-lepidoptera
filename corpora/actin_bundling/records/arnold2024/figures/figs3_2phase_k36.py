"""Calibrated reader for SI Fig. S3 (left panel, 2-phase bilayer; PDF p. 13 / SI-4):
Actin DDM image structure function D(k, t) vs. time t [s] at k = 3.6 um^-1 (turquoise/teal, drawn on top).
"""

from __future__ import annotations

import numpy as np

from boxgraph import digitize as D

CURVE_INDEX = 3  # 0: k=0.9, 1: k=1.8, 2: k=2.7, 3: k=3.6 um^-1
Y_TICKS_TOP_DOWN = [3.0e10, 2.5e10, 2.0e10, 1.5e10, 1.0e10, 0.5e10]


def _group(indices: np.ndarray, gap: int = 2) -> list[list[int]]:
    out: list[list[int]] = []
    for x in indices:
        xi = int(x)
        if out and xi <= out[-1][-1] + gap:
            out[-1].append(xi)
        else:
            out.append([xi])
    return out


def _calibrate_and_colors(image: np.ndarray):
    H, W = image.shape[:2]
    r, g, b = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    dark = (mx < 140) & ((mx - mn) < 30)

    hg_top = _group(np.where(dark[:400, :].sum(axis=1) > 0.5 * W)[0], gap=4)[0]
    vg = _group(np.where(dark[:, 200:].sum(axis=0) > 0.5 * H)[0] + 200, gap=4)
    r_top_mid = (hg_top[0] + hg_top[-1]) / 2.0 + 0.5
    r_top_in = hg_top[-1]
    c_left_mid = (vg[0][0] + vg[0][-1]) / 2.0 + 0.5
    c_left_in = vg[0][-1]
    c_right_mid = (vg[-1][0] + vg[-1][-1]) / 2.0 + 0.5
    c_right_in = vg[-1][0]

    tx_raw = []
    for c in range(c_left_in + 20, c_right_in - 20):
        n = 0
        while r_top_in + 1 + n < H and dark[r_top_in + 1 + n, c]:
            n += 1
        if 20 <= n <= 45:
            tx_raw.append(c)
    tx_pos = sorted((grp[0] + grp[-1]) / 2.0 + 0.5 for grp in _group(np.array(tx_raw), gap=1) if len(grp) <= 10)
    if len(tx_pos) != 2:
        raise D.CalibrationError(f"Fig. S3 left: expected 2 interior major x-ticks, got {tx_pos}")
    xa = D.Axis("log10", [(c_left_mid, 10.0), (tx_pos[0], 100.0), (tx_pos[1], 1000.0), (c_right_mid, 10000.0)])

    ry_raw = []
    for r_idx in range(r_top_in + 40, min(H - 10, r_top_in + 1800)):
        n = 0
        while c_right_in - 1 - n >= 0 and dark[r_idx, c_right_in - 1 - n]:
            n += 1
        if 15 <= n <= 45:
            ry_raw.append(r_idx)
    ry_pos = sorted((grp[0] + grp[-1]) / 2.0 + 0.5 for grp in _group(np.array(ry_raw), gap=1) if len(grp) <= 10)
    if len(ry_pos) != len(Y_TICKS_TOP_DOWN) - 1:
        raise D.CalibrationError(f"Fig. S3 left: expected {len(Y_TICKS_TOP_DOWN) - 1} right y-ticks, got {ry_pos}")
    ya = D.Axis("linear", list(zip([r_top_mid] + ry_pos, Y_TICKS_TOP_DOWN)))

    c_leg0 = int(c_left_mid) + 120
    c_leg1 = int(c_left_mid) + 340
    leg = image[260:900, c_leg0:c_leg1]
    lsat = leg.max(axis=2).astype(np.int16) - leg.min(axis=2).astype(np.int16)
    cg = [grp for grp in _group(np.where((lsat > 40).any(axis=1))[0], gap=10) if len(grp) >= 15]
    if len(cg) != 4:
        raise D.CalibrationError(f"Fig. S3 left: expected 4 legend color entries, got {len(cg)}")
    colors = []
    for grp in cg:
        band = leg[grp[0] : grp[-1] + 1]
        bsat = band.max(axis=2).astype(np.int16) - band.min(axis=2).astype(np.int16)
        top = bsat >= np.percentile(bsat[bsat > 40], 80)
        colors.append(np.median(band[top], axis=0).astype(np.float32))

    return xa, ya, colors, (r_top_in, int(round(ya.to_pixel(0.0))), c_left_in, c_right_in)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, colors, (r_top_in, r_zero, c_left_in, c_right_in) = _calibrate_and_colors(image)
    H, W = image.shape[:2]

    rgb = image.astype(np.float32)
    target = colors[CURVE_INDEX]
    d_col = np.sqrt(((rgb - target) ** 2).sum(axis=2))
    mask = d_col < 18.0
    mask[: r_top_in + 2, :] = False
    mask[min(H, r_zero + 30) :, :] = False
    mask[:, : c_left_in + 2] = False
    mask[:, c_right_in - 2 :] = False
    mask[220:920, c_left_in + 50 : c_left_in + 1200] = False

    base_err_px = 3.0 + max(ya.residual_px, 0.5)
    out = []
    for xv in at:
        pc = float(xa.to_pixel(xv))
        c0 = int(round(pc))
        hw = 18
        c_lo, c_hi = max(c_left_in + 2, c0 - hw), min(c_right_in - 2, c0 + hw + 1)
        rr = np.where(mask[:, c_lo:c_hi])[0]
        if len(rr) < 20:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no k=3.6 pixels near t={xv:g}"})
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
            "note": f"2-phase k=3.6 um^-1 (10-90 pct + cal); cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px",
        })
    return out
