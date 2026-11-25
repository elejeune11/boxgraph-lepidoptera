"""Calibrated reader for Fig. 4C (PDF p. 6 / printed p. 26575):
Mean tracer filament position variance <|Delta x|^2> - <|Delta x|>^2 [um^2] vs. time t [min],
with SEM error bars (n = 212 in caption / 212-254 in text).

Note on y-axis typesetting: the five MATLAB y-axis ticks sit at exact 452.5 px intervals
(0, 1, 2, 3, 4 x 10^-3 um^2, residual 0.00 px), whereas the manually re-typeset digit labels
'3', '2', '1' drift downward from their ticks. Calibration uses the exact axis tick marks.
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
    black = (mx < 120) & ((mx - mn) < 25)

    # Bottom x-axis spine and left y-axis spine
    bot_g = _group(np.where(black.sum(axis=1) > 0.5 * W)[0], gap=4)[-1]
    left_g = _group(np.where(black.sum(axis=0) > 0.5 * H)[0], gap=4)[0]
    r_x0, r_x1 = bot_g[0], bot_g[-1]
    c_y0, c_y1 = left_g[0], left_g[-1]
    r_zero = (r_x0 + r_x1) / 2.0 + 0.5
    c_zero = (c_y0 + c_y1) / 2.0 + 0.5

    # X-axis ticks pointing upward from bottom spine: 5, 10, 15 min
    xt_raw = []
    for c in range(c_y1 + 10, W - 5):
        n = 0
        while r_x0 - 1 - n >= 0 and black[r_x0 - 1 - n, c]:
            n += 1
        if 12 <= n <= 60:
            xt_raw.append(c)
    xt_groups = [g for g in _group(np.array(xt_raw), gap=1) if 4 <= len(g) <= 20]
    xt_pos = sorted((g[0] + g[-1]) / 2.0 + 0.5 for g in xt_groups)
    # Match to 0 (left spine) + 5, 10, 15
    x_ticks = D.regular_ticks([c_zero] + xt_pos, tol=3.0)
    if len(x_ticks) != 4:
        raise D.CalibrationError(f"Fig. 4C: expected 4 x-ticks (0,5,10,15), got {x_ticks}")
    xa = D.Axis("linear", list(zip(x_ticks, [0.0, 5.0, 10.0, 15.0])))

    # Y-axis ticks pointing rightward from left spine (length ~24 px, thickness ~9-10 px;
    # note the t=0 error bar caps also protrude to the left of the left spine, whereas
    # genuine y-axis ticks protrude ONLY to the right of the left spine!)
    yt_raw = []
    for r_idx in range(50, r_x0 - 10):
        # Genuine ticks do NOT have black pixels immediately left of c_y0
        if black[r_idx, max(0, c_y0 - 12):c_y0 - 2].any():
            continue
        n = 0
        while c_y1 + 1 + n < W and black[r_idx, c_y1 + 1 + n]:
            n += 1
        if 15 <= n <= 40:
            yt_raw.append(r_idx)
    yt_groups = [g for g in _group(np.array(yt_raw), gap=1) if 5 <= len(g) <= 16]
    yt_pos = sorted((g[0] + g[-1]) / 2.0 + 0.5 for g in yt_groups)
    y_ticks = D.regular_ticks(yt_pos + [r_zero], tol=3.0)
    if len(y_ticks) != 5:
        raise D.CalibrationError(f"Fig. 4C: expected 5 y-ticks (4,3,2,1,0), got {y_ticks}")
    # Topmost row is 4e-3, bottommost is 0
    ya = D.Axis("linear", list(zip(y_ticks, [4e-3, 3e-3, 2e-3, 1e-3, 0.0])))

    return xa, ya, black, (c_y0, c_y1, r_x0, r_x1)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, black, (c_y0, c_y1, r_x0, r_x1) = _calibrate(image)
    H, W = image.shape[:2]

    # Restrict to plot area (exclude the 'Tracer trajectory' text/inset in the upper middle).
    # Note the inset micrograph is colored, while the title 'Tracer trajectory' sits at rows < 360
    # and cols between x=2.5 and x=12 min.
    c_inset_l = int(round(xa.to_pixel(2.5)))
    c_inset_r = int(round(xa.to_pixel(11.0)))
    r_inset_b = int(round(ya.to_pixel(1.55e-3)))

    clean = black.copy()
    clean[:150, :] = False
    clean[r_x0 - 2:, :] = False
    clean[:, :c_y0 - 40] = False
    clean[150:r_inset_b, c_inset_l:c_inset_r] = False

    # Erase left spine EXCEPT where wider horizontal features (t=0 marker & error bar caps) stick out
    # Actually, morphological opening with a 19x19 circular disk removes ALL ~9-10 px spines,
    # connecting lines, error-bar stems, and error-bar caps, leaving ONLY the ~32 px circular markers!
    yy, xx = np.ogrid[-9:10, -9:10]
    disk = (xx * xx + yy * yy) <= 9 * 9
    markers_mask = ndimage.binary_opening(clean, structure=disk)

    lab, n_lab = ndimage.label(markers_mask)
    pts = []
    for i in range(1, n_lab + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 300:  # a 32-px circle has area ~800 px; after opening still ~600-800 px
            continue
        r_mid = float(rr.mean()) + 0.5
        c_mid = float(cc.mean()) + 0.5
        if c_mid < c_y0 - 15 or c_mid > xa.to_pixel(15.2):
            continue
        pts.append((c_mid, r_mid, int(rr.min()), int(rr.max())))
    pts.sort(key=lambda p: p[0])

    # For each detected marker, locate its upper and lower horizontal error bar caps
    # (caps are ~46 px wide and ~8-9 px tall, centered on c_mid)
    y_misfit = max(ya.residual_px, 0.5)
    cap_err_px = 1.5 + y_misfit
    pt_err_px = 1.5 + y_misfit

    out = []
    for c_mid, r_mid, r_top_m, r_bot_m in pts:
        c0 = int(round(c_mid))
        # Search vertically along a band centered at c0 of half-width 18 px (width 37 px < 46 px cap)
        c_lo, c_hi = max(0, c0 - 18), min(W, c0 + 19)
        col_band = clean[:, c_lo:c_hi]
        # A horizontal cap fills >85% of this 37-px wide band!
        cap_rows = col_band.mean(axis=1) > 0.85
        # Find cap above r_top_m - 5 and cap below r_bot_m + 5
        above_idx = np.where(cap_rows[: max(0, r_top_m - 4)])[0]
        below_idx = np.where(cap_rows[r_bot_m + 5 : r_x0 - 4])[0] + (r_bot_m + 5)

        above_groups = _group(above_idx, gap=2)
        below_groups = _group(below_idx, gap=2)
        if not above_groups or not below_groups:
            raise RuntimeError(f"Fig. 4C: missing error bar cap for marker at col={c_mid:.1f}, row={r_mid:.1f}")
        # Take the cap closest to the marker on each side (for t=0, ticks don't span left of spine,
        # so they fail the 37-px centered band check anyway!)
        top_cap = above_groups[-1]
        bot_cap = below_groups[0]
        r_bar_top = (top_cap[0] + top_cap[-1]) / 2.0 + 0.5
        r_bar_bot = (bot_cap[0] + bot_cap[-1]) / 2.0 + 0.5

        x_val = float(xa.to_value(c_mid))
        y_val = float(ya.to_value(r_mid))
        y_lo = float(min(ya.to_value(r_mid - pt_err_px), ya.to_value(r_mid + pt_err_px)))
        y_hi = float(max(ya.to_value(r_mid - pt_err_px), ya.to_value(r_mid + pt_err_px)))

        b_hi = float(ya.to_value(r_bar_top))
        b_lo = float(ya.to_value(r_bar_bot))
        b_hi_out = float(ya.to_value(r_bar_top - cap_err_px))
        b_lo_out = float(ya.to_value(r_bar_bot + cap_err_px))

        out.append({
            "x": x_val,
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "bar_lo": min(b_lo, b_hi),
            "bar_hi": max(b_lo, b_hi),
            "bar_lo_outer": min(b_lo_out, b_hi_out),
            "bar_hi_outer": max(b_lo_out, b_hi_out),
            "px": [c_mid, r_mid],
            "px_lo": float(r_bar_top),
            "px_hi": float(r_bar_bot),
            "note": (
                f"SEM error bar [{min(b_lo, b_hi):.4g}, {max(b_lo, b_hi):.4g}] um^2; "
                f"cal misfit x={xa.residual_px:.2f} px, y={ya.residual_px:.2f} px"
            ),
        })
    return out
