"""Shared calibrated reader for Extended Data Fig. 5c wetting tension profiles gamma(e).

The vertical axis linearly spans gamma_AVR (= 1e-8 pN/nm, bottom dashed line and tick 0)
to gamma_AVD (= 1e4 pN/nm, top dashed line and tick 4) with 5 evenly spaced ticks
corresponding to normalized fractions 0.00, 0.25, 0.50, 0.75, 1.00.
"""

from __future__ import annotations

import numpy as np
import _calib as C

GAMMA_AVR = 1.0e-8
GAMMA_AVD = 1.0e4

# Legend order top-to-bottom: f(e) = e-1, 8(e-1), 16(e-1), 32(e-1), 64(e-1)
CURVE_NAMES = ("ke1", "ke8", "ke16", "ke32", "ke64")


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 10, int(vr0[2]) + 10), "x", side=+1, band=(4, 14), frac=0.4)
    yt = C.find_outward_ticks(dark, vl0, (int(ht[1]) - 10, int(hb[2]) + 10), "y", side=-1, band=(4, 14), frac=0.4)
    xt = C.match_ticks_regular(xt, 5, tol=4.0)
    yt = C.match_ticks_regular(yt, 5, tol=3.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 0.5, 1.0, 1.5, 2.0])), max_residual_px=2.0)
    # 5 evenly spaced y-ticks from gamma_AVD (top) down to gamma_AVR (bottom)
    y_vals = [GAMMA_AVR + frac * (GAMMA_AVD - GAMMA_AVR) for frac in (1.0, 0.75, 0.50, 0.25, 0.0)]
    ya = C.verify_axis("linear", list(zip(yt, y_vals)), max_residual_px=1.5)

    # Sample the 5 curve RGB colors directly from the legend box inside the image
    c_leg0 = int(xa.to_pixel(1.06))
    c_leg1 = int(xa.to_pixel(1.38))
    r_leg0 = int(ya.to_pixel(0.55 * GAMMA_AVD))
    r_leg1 = int(ya.to_pixel(0.05 * GAMMA_AVD))
    leg = image[r_leg0:r_leg1, c_leg0:c_leg1].astype(float)
    sat = np.maximum(np.maximum(leg[:, :, 0], leg[:, :, 1]), leg[:, :, 2]) - np.minimum(
        np.minimum(leg[:, :, 0], leg[:, :, 1]), leg[:, :, 2]
    )
    bands = C.clusters(np.where((sat > 40).mean(axis=1) > 0.35)[0] + r_leg0, gap=6)
    if len(bands) != 5:
        raise RuntimeError(f"Expected 5 legend bands in ED Fig. 5c, found {len(bands)}")
    colors = {}
    for name, (_, b0, b1) in zip(CURVE_NAMES, bands):
        patch = image[b0 : b1 + 1, c_leg0 + 5 : c_leg1 - 5].astype(float)
        psat = np.maximum(np.maximum(patch[:, :, 0], patch[:, :, 1]), patch[:, :, 2]) - np.minimum(
            np.minimum(patch[:, :, 0], patch[:, :, 1]), patch[:, :, 2]
        )
        colors[name] = np.median(patch[psat >= 0.85 * float(psat.max())], axis=0)

    row_min = int(yt[0]) - 18
    row_max = int(yt[-1]) + 18
    leg_box = (int(ya.to_pixel(0.50 * GAMMA_AVD)), int(ya.to_pixel(0.03 * GAMMA_AVD)), int(xa.to_pixel(1.06)), int(vr0[1]) - 5)
    return xa, ya, colors, row_min, row_max, leg_box


def read_curve(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, colors, row_min, row_max, leg_box = calibrate(image)
    rgb = image.astype(float)
    dists = {name: np.linalg.norm(rgb - col[None, None, :], axis=2) for name, col in colors.items()}
    d_target = dists[which]
    d_other = np.minimum.reduce([d for name, d in dists.items() if name != which])
    # Require both close to target color AND strictly closer than any other curve color
    mask = (d_target < 55.0) & (d_target + 15.0 < d_other)

    # Mask legend box so legend lines are never mistaken for curves
    lr0, lr1, lc0, lc1 = leg_box
    mask[lr0:lr1, lc0:lc1] = False

    # Mask vertical red dashed interface line at e = 1.0
    c_int = int(round(xa.to_pixel(1.0)))
    mask[:, max(0, c_int - 7) : c_int + 8] = False

    # Mask top blue dashed reference line gamma_AVD for ke1 (ke1 stays below 0.9*gamma_AVD on [0, 2])
    r_top = int(round(ya.to_pixel(GAMMA_AVD)))
    r_bot = int(round(ya.to_pixel(GAMMA_AVR)))
    if which == "ke1":
        mask[max(0, r_top - 12) : r_top + 13, :] = False

    out = []
    for x in at:
        res = C.sample_solid_curve_y(mask, xa, ya, x, row_min, row_max, search_dx=4, note_prefix=f"ED5c {which}: ")
        if (not res.get("found", res.get("y") is not None)) and which != "ke1" and abs(x - 1.0) > 0.05:
            # Flat plateau regions where later curves (purple ke64) paint over earlier ones
            c_f = float(xa.to_pixel(x))
            c = int(round(c_f))
            r_ref = r_bot if x < 1.0 else r_top
            band = image[r_ref - 12 : r_ref + 13, max(0, c - 4) : c + 5].astype(float)
            nonwhite = (band.mean(axis=2) < 230).any(axis=1)
            rs = C.runs_1d(nonwhite)
            if rs:
                r0, r1 = max(rs, key=lambda t: t[1] - t[0])
                r_mid = (r_ref - 12) + 0.5 * (r0 + r1) + 0.5
                half_px = max(0.5 * (r1 - r0 + 1), 2.5) + 1.0 + ya.residual_px
                yv = float(ya.to_value(r_mid))
                ya_, yb_ = float(ya.to_value(r_mid - half_px)), float(ya.to_value(r_mid + half_px))
                lo, hi = min(ya_, yb_), max(ya_, yb_)
                res = {
                    "x": float(x),
                    "y": yv,
                    "lo": lo,
                    "hi": hi,
                    "px": [c_f, float(r_mid)],
                    "px_lo": float(r_mid - half_px),
                    "px_hi": float(r_mid + half_px),
                    "note": f"ED5c {which}: overlapping plateau, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
        out.append(res)
    return out
