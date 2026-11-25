"""Calibrated reader for Fig. 3B (PDF p. 5 / printed p. 26574):
Actin intermediate scattering function F(<k>, t) vs. time t [s] on single-phase (Ld-only) membrane
(blue/cyan open triangles, t = 100 to 6500 s).
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
    dark = (mx < 140) & ((mx - mn) < 35)

    hg = _group(np.where(dark.sum(axis=1) > 0.3 * W)[0], gap=4)
    r_top_mid = (hg[0][0] + hg[0][-1]) / 2.0 + 0.5
    r_bot_mid = (hg[-1][0] + hg[-1][-1]) / 2.0 + 0.5
    r_top_in = hg[0][-1]
    r_bot_in = hg[-1][0]

    vg = _group(np.where(dark[:, 500:].sum(axis=0) > 0.3 * H)[0] + 500, gap=4)
    c_left_mid = (vg[0][0] + vg[0][-1]) / 2.0 + 0.5
    c_right_mid = (vg[-1][0] + vg[-1][-1]) / 2.0 + 0.5
    c_left_in = vg[0][-1]
    c_right_in = vg[-1][0]

    c_mid_est = int(round((c_left_mid + c_right_mid) / 2.0))
    xt_mid = []
    for c in range(c_mid_est - 80, c_mid_est + 80):
        n = 0
        while r_bot_in - 1 - n >= 0 and dark[r_bot_in - 1 - n, c]:
            n += 1
        if 15 <= n <= 35:
            xt_mid.append(c)
    mid_groups = [grp for grp in _group(np.array(xt_mid), gap=1) if len(grp) <= 20]
    if not mid_groups:
        raise D.CalibrationError("Fig. 3B: could not find major x-tick at 10^3")
    c_1000 = (mid_groups[0][0] + mid_groups[0][-1]) / 2.0 + 0.5
    xa = D.Axis("log10", [(c_left_mid, 100.0), (c_1000, 1000.0), (c_right_mid, 10000.0)])

    yt_raw = []
    for r_idx in range(r_top_in + 40, r_bot_in - 40):
        n = 0
        while c_right_in - 1 - n >= 0 and dark[r_idx, c_right_in - 1 - n]:
            n += 1
        if 12 <= n <= 40:
            yt_raw.append(r_idx)
    yt_groups = [grp for grp in _group(np.array(yt_raw), gap=1) if len(grp) <= 20]
    yt_pos = sorted((grp[0] + grp[-1]) / 2.0 + 0.5 for grp in yt_groups)
    y_ticks = D.regular_ticks([r_top_mid] + yt_pos + [r_bot_mid], tol=6.0)
    if len(y_ticks) != 7:
        raise D.CalibrationError(f"Fig. 3B: expected 7 y-ticks (1.0..-0.2), got {y_ticks}")
    ya = D.Axis("linear", list(zip(y_ticks, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2])))

    return xa, ya, (r_top_in, r_bot_in, c_left_in, c_right_in)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, (r_top_in, r_bot_in, c_left_in, c_right_in) = _calibrate(image)
    H, W = image.shape[:2]

    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    blue = (b - r > 50) & (g - r > 20) & (b > 100)

    mask = np.zeros((H, W), dtype=bool)
    mask[max(0, r_top_in - 45) : r_bot_in - 2, max(0, c_left_in - 20) : c_right_in] = blue[
        max(0, r_top_in - 45) : r_bot_in - 2, max(0, c_left_in - 20) : c_right_in
    ]
    # Exclude the blue 'Ld only' text label (t > 3200, F > 0.72)
    c_lbl0 = int(round(xa.to_pixel(3200.0)))
    r_lbl1 = int(round(ya.to_pixel(0.72)))
    mask[r_top_in:r_lbl1, c_lbl0:c_right_in] = False

    base_err_px = 5.0 + max(ya.residual_px, 1.0)

    out = []
    for xv in at:
        pc = float(xa.to_pixel(xv))
        c0 = int(round(pc))
        hw = 18 if xv < 4500.0 else 28
        c_lo, c_hi = max(0, c0 - hw), min(W, c0 + hw + 1)
        rr = np.where(mask[:, c_lo:c_hi])[0]
        if len(rr) < 15:
            c_lo, c_hi = max(0, c0 - 35), min(W, c0 + 36)
            rr = np.where(mask[:, c_lo:c_hi])[0]
        if len(rr) < 15:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no blue pixels near t={xv:g}"})
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
            "note": f"blue triangles (10-90 pct scatter + cal); cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px",
        })
    return out
