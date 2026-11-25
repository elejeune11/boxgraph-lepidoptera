"""Calibrated reader for Fig. 4B (bottom panel; PDF p. 6 / printed p. 26575):
Domain excess perimeter L_ex vs. time t [min] during cholesterol insertion (t = 0 to 11 min),
with error bars (type and n unstated in paper).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

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
    black = (mx < 130) & ((mx - mn) < 35)

    hg = _group(np.where(black.sum(axis=1) > 0.5 * W)[0], gap=6)
    r_top_x0, r_top_x1 = hg[0][0], hg[0][-1]
    r_bot_x0, r_bot_x1 = hg[1][0], hg[1][-1]

    # Top-left spine aligns with x=0 ticks; bottom-left spine carries the bottom y ticks
    tl_g = _group(np.where(black[:r_top_x0, 600:1000].sum(axis=0) > 0.5 * r_top_x0)[0] + 600, gap=4)[0]
    bl_g = _group(
        np.where(black[r_top_x1 + 50 : r_bot_x0, 600:1000].sum(axis=0) > 0.5 * (r_bot_x0 - r_top_x1 - 50))[0] + 600,
        gap=4,
    )[0]
    c_zero = (tl_g[0] + tl_g[-1]) / 2.0 + 0.5
    c_bl0, c_bl1 = bl_g[0], bl_g[-1]

    # Bottom x-axis ticks: 0, 5, 10, 15, 20 min
    xt_raw = []
    for c in range(c_bl1 + 20, W - 10):
        n = 0
        while r_bot_x0 - 1 - n >= 0 and black[r_bot_x0 - 1 - n, c]:
            n += 1
        if 14 <= n <= 60:
            xt_raw.append(c)
    xt_groups = [g for g in _group(np.array(xt_raw), gap=1) if 8 <= len(g) <= 30]
    xt_pos = sorted((g[0] + g[-1]) / 2.0 + 0.5 for g in xt_groups)
    x_ticks = D.regular_ticks([c_zero] + xt_pos, tol=4.0)
    if len(x_ticks) != 5:
        raise D.CalibrationError(f"Fig. 4B bottom L_ex: expected 5 x-ticks, got {x_ticks}")
    xa = D.Axis("linear", list(zip(x_ticks, [0.0, 5.0, 10.0, 15.0, 20.0])))

    # Bottom panel left y-axis ticks (2.0, 1.5, 1.0, 0.5; pointing rightward from bottom-left spine)
    yt_raw = []
    for r_idx in range(r_top_x1 + 50, r_bot_x0 - 10):
        # Skip the t=0 marker/cap which also protrudes left of c_bl0
        if black[r_idx, max(0, c_bl0 - 20) : c_bl0 - 2].any():
            continue
        n = 0
        while c_bl1 + 1 + n < W and black[r_idx, c_bl1 + 1 + n]:
            n += 1
        if 25 <= n <= 70:
            yt_raw.append(r_idx)
    yt_groups = [g for g in _group(np.array(yt_raw), gap=1) if 10 <= len(g) <= 30]
    yt_pos = sorted((g[0] + g[-1]) / 2.0 + 0.5 for g in yt_groups)
    y_ticks = D.regular_ticks(yt_pos, tol=4.0)
    if len(y_ticks) != 4:
        raise D.CalibrationError(f"Fig. 4B bottom L_ex: expected 4 y-ticks (2.0,1.5,1.0,0.5), got {y_ticks}")
    ya = D.Axis("linear", list(zip(y_ticks, [2.0, 1.5, 1.0, 0.5])))

    return xa, ya, black, (c_bl0, c_bl1, r_top_x1, r_bot_x0)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, black, (c_bl0, c_bl1, r_top_x1, r_bot_x0) = _calibrate(image)
    H, W = image.shape[:2]

    # Restrict to bottom panel left of dashed line (t < 12.2 min), excluding the '2-phase' text
    # which sits in the lower right (t > 5.5 min, L_ex < 0.75)
    c_right_lim = int(round(xa.to_pixel(12.2)))
    clean = np.zeros((H, W), dtype=bool)
    clean[r_top_x1 + 100 : r_bot_x0 - 10, max(0, c_bl0 - 70) : c_right_lim] = black[
        r_top_x1 + 100 : r_bot_x0 - 10, max(0, c_bl0 - 70) : c_right_lim
    ]
    # Mask out the '2-phase' label region (x > 5.5 min and y < 0.85)
    c_lbl = int(round(xa.to_pixel(5.5)))
    r_lbl = int(round(ya.to_pixel(0.85)))
    clean[r_lbl:r_bot_x0, c_lbl:c_right_lim] = False

    # Isolate the 9 oval markers with fast separable 1D box erosion (29 x 55 px) on the panel ROI
    r0_roi, r1_roi = r_top_x1 + 100, r_bot_x0 - 10
    c0_roi, c1_roi = max(0, c_bl0 - 70), c_right_lim
    roi = clean[r0_roi:r1_roi, c0_roi:c1_roi]
    er = ndimage.minimum_filter1d(roi, size=29, axis=0)
    er = ndimage.minimum_filter1d(er, size=55, axis=1)

    lab, n_lab = ndimage.label(er)
    pts = []
    for i in range(1, n_lab + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 100:
            continue
        r_mid = float(rr.mean()) + 0.5 + r0_roi
        c_mid = float(cc.mean()) + 0.5 + c0_roi
        pts.append((c_mid, r_mid))
    pts.sort(key=lambda p: p[0])

    if len(pts) != 9:
        raise RuntimeError(f"Fig. 4B L_ex: expected 9 markers, found {len(pts)}")

    # Error bar caps are ~90 px wide and ~20 px tall, centered at each marker's x column.
    # For t = 0..5 min the caps sit directly against the oval marker, while for t = 6, 8, 11 min
    # white space separates cap and marker. In all 9 cases, requiring >85% fill across a 75-px wide
    # horizontal band centered at c_mid cleanly selects the marker + both caps, whose outermost
    # top and bottom edges lie 9.5 px outside the top and bottom cap centerlines.
    y_misfit = max(ya.residual_px, 0.5)
    cap_err_px = 2.0 + y_misfit
    pt_err_px = 2.0 + y_misfit

    out = []
    for c_mid, r_mid in pts:
        c0 = int(round(c_mid))
        band = clean[:, max(0, c0 - 37) : min(W, c0 + 38)]
        # Search within ±200 px vertically of the marker center (max error bar half-height is ~70 px)
        r_lo_s = max(r_top_x1 + 100, int(round(r_mid)) - 200)
        r_hi_s = min(r_bot_x0 - 10, int(round(r_mid)) + 200)
        solid = band[r_lo_s:r_hi_s, :].mean(axis=1) > 0.85
        rows = np.where(solid)[0] + r_lo_s
        if len(rows) < 20:
            raise RuntimeError(f"Fig. 4B L_ex: failed to find error bar caps at col={c_mid:.1f}")
        r_bar_top = float(rows[0]) + 9.5
        r_bar_bot = float(rows[-1]) - 9.5
        # Marker center from the symmetric error-bar midpoint (avoids connecting-line pull);
        # check it agrees closely with r_mid
        r_ctr = (r_bar_top + r_bar_bot) / 2.0
        if abs(r_ctr - r_mid) > 8.0:
            r_ctr = r_mid

        x_val = float(xa.to_value(c_mid))
        y_val = float(ya.to_value(r_ctr))
        a, b = ya.to_value(r_ctr - pt_err_px), ya.to_value(r_ctr + pt_err_px)
        b_hi = float(ya.to_value(r_bar_top))
        b_lo = float(ya.to_value(r_bar_bot))
        b_hi_out = float(ya.to_value(r_bar_top - cap_err_px))
        b_lo_out = float(ya.to_value(r_bar_bot + cap_err_px))

        out.append({
            "x": x_val,
            "y": y_val,
            "lo": float(min(a, b)),
            "hi": float(max(a, b)),
            "bar_lo": min(b_lo, b_hi),
            "bar_hi": max(b_lo, b_hi),
            "bar_lo_outer": min(b_lo_out, b_hi_out),
            "bar_hi_outer": max(b_lo_out, b_hi_out),
            "px": [c_mid, r_ctr],
            "px_lo": float(r_bar_top),
            "px_hi": float(r_bar_bot),
            "note": (
                f"unstated error bar [{min(b_lo, b_hi):.4g}, {max(b_lo, b_hi):.4g}]; "
                f"cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px"
            ),
        })
    return out
