"""Calibrated reader for Fig. 3C (PDF p. 5 / printed p. 26574):
Mean wavelength a = 2*pi/<k> [um] of Lo lipid domains vs. time t [s] (green squares, t = 0 to 4000 s).
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
    dark = (mx < 130) & ((mx - mn) < 35)

    # Top and bottom horizontal spines; left and right vertical spines
    hg = _group(np.where(dark.sum(axis=1) > 0.4 * W)[0], gap=4)
    vg = _group(np.where(dark.sum(axis=0) > 0.4 * H)[0], gap=4)
    r_top_mid = (hg[0][0] + hg[0][-1]) / 2.0 + 0.5
    r_bot_mid = (hg[-1][0] + hg[-1][-1]) / 2.0 + 0.5
    r_top_in = hg[0][-1]
    r_bot_in = hg[-1][0]
    c_left_in = vg[0][-1]
    c_right_in = vg[-1][0]
    c_right_mid = (vg[-1][0] + vg[-1][-1]) / 2.0 + 0.5

    # X-axis ticks along bottom spine (1000, 2000, 3000) + right spine (4000)
    # Note: left spine sits ~5 px left of t=0; the 3 interior ticks + right spine have 0.00 px misfit.
    xt_raw = []
    for c in range(c_left_in + 10, c_right_in - 10):
        n = 0
        while r_bot_in - 1 - n >= 0 and dark[r_bot_in - 1 - n, c]:
            n += 1
        if 12 <= n <= 50:
            xt_raw.append(c)
    xt_groups = [grp for grp in _group(np.array(xt_raw), gap=1) if len(grp) <= 20]
    xt_pos = sorted((grp[0] + grp[-1]) / 2.0 + 0.5 for grp in xt_groups)
    x_ticks = D.regular_ticks(xt_pos + [c_right_mid], tol=4.0)
    if len(x_ticks) != 4:
        raise D.CalibrationError(f"Fig. 3C: expected 4 x-ticks (1000,2000,3000,4000), got {x_ticks}")
    xa = D.Axis("linear", list(zip(x_ticks, [1000.0, 2000.0, 3000.0, 4000.0])))

    # Y-axis ticks along right spine (unobstructed: 3.5, 3.0) + top (4.0) and bottom (2.5) spines
    yt_raw = []
    for r_idx in range(r_top_in + 10, r_bot_in - 10):
        n = 0
        while c_right_in - 1 - n >= 0 and dark[r_idx, c_right_in - 1 - n]:
            n += 1
        if 12 <= n <= 50:
            yt_raw.append(r_idx)
    yt_groups = [grp for grp in _group(np.array(yt_raw), gap=1) if len(grp) <= 20]
    yt_pos = sorted((grp[0] + grp[-1]) / 2.0 + 0.5 for grp in yt_groups)
    y_ticks = D.regular_ticks([r_top_mid] + yt_pos + [r_bot_mid], tol=4.0)
    if len(y_ticks) != 4:
        raise D.CalibrationError(f"Fig. 3C: expected 4 y-ticks (4.0,3.5,3.0,2.5), got {y_ticks}")
    ya = D.Axis("linear", list(zip(y_ticks, [4.0, 3.5, 3.0, 2.5])))

    return xa, ya, (r_top_in, r_bot_in, c_left_in, c_right_in)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, (r_top_in, r_bot_in, c_left_in, c_right_in) = _calibrate(image)
    H, W = image.shape[:2]

    r = image[:, :, 0].astype(np.int16)
    g = image[:, :, 1].astype(np.int16)
    b = image[:, :, 2].astype(np.int16)
    green = (g - r > 40) & (g - b > 30) & (g > 70)
    purple = (r - g > 30) & (b - g > 35)

    # Restrict to plot interior (allow slight overhang at t=0 and t=4000 where markers straddle spines)
    mask = np.zeros((H, W), dtype=bool)
    mask[r_top_in + 2 : r_bot_in - 2, max(0, c_left_in - 30) : min(W, c_right_in + 35)] = green[
        r_top_in + 2 : r_bot_in - 2, max(0, c_left_in - 30) : min(W, c_right_in + 35)
    ]
    # Exclude the green text label 'Domain a' (t in [700, 2100], y > 3.55; green curve stays below 3.48 for t <= 2100)
    c_lbl0 = int(round(xa.to_pixel(700.0)))
    c_lbl1 = int(round(xa.to_pixel(2100.0)))
    r_lbl1 = int(round(ya.to_pixel(3.55)))
    mask[r_top_in:r_lbl1, c_lbl0:c_lbl1] = False

    y_misfit = max(ya.residual_px, 1.0)
    # 1500 dpi is 5x the 300 dpi embedded plot raster -> add >= 5 px (1 native raster pixel) + calibration
    base_err_px = 5.0 + y_misfit

    out = []
    for xv in at:
        pc = float(xa.to_pixel(xv))
        c0 = int(round(pc))
        hw = 15  # +-25 s window (~half a marker width)
        c_lo, c_hi = max(0, c0 - hw), min(W, c0 + hw + 1)
        rr = np.where(mask[:, c_lo:c_hi])[0]
        if len(rr) < 20:
            out.append({"x": float(xv), "y": None, "found": False, "note": f"no green pixels near t={xv:g}"})
            continue
        r_med = float(np.median(rr)) + 0.5
        r_p10 = float(np.percentile(rr, 10)) + 0.5
        r_p90 = float(np.percentile(rr, 90)) + 0.5

        # Check if purple actin curve overlaps in the same vertical band (crossing around t = 800..1150 s);
        # if so, widen the vertical interval to account for partial occlusion of the green square markers.
        purp_rr = np.where(purple[r_top_in + 2 : r_bot_in - 2, c_lo:c_hi])[0] + (r_top_in + 2)
        occluded = False
        if len(purp_rr) > 20 and abs(float(np.median(purp_rr)) - r_med) < 75.0:
            occluded = True
            pad = 22.0 + base_err_px
        else:
            pad = base_err_px

        r_hi_edge = min(r_p10, r_med) - pad
        r_lo_edge = max(r_p90, r_med) + pad
        y_val = float(ya.to_value(r_med))
        a, b = ya.to_value(r_hi_edge), ya.to_value(r_lo_edge)
        note = (
            f"green squares (partially overlapped by purple actin curve); cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px"
            if occluded
            else f"green squares; cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px"
        )
        out.append({
            "x": float(xv),
            "y": y_val,
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "px": [pc, r_med],
            "px_lo": float(r_hi_edge),
            "px_hi": float(r_lo_edge),
            "note": note,
        })
    return out
