"""Calibrated reader for Fig. 2c (p. 4 / 450): fitted 2*k_on/D (mislabeled k_on/D on y-axis), 5 red squares + ±20% bars."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import NATIVE_QUANT_PX, _group_runs, calibrate_fig2c, find_hline_row_red, value_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2c(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0

    # Detect the 5 orange-filled red squares on 1/4-res:
    # solid red 8x8 block with substantial orange interior fill (G - B > 15).
    sub = image[::4, ::4].astype(float)
    sR, sG, sB = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    sint = np.zeros(sR.shape, dtype=bool)
    sint[int(f["tr1"]) // 4 + 3 : int(f["br0"]) // 4 - 3, int(f["lc1"]) // 4 + 3 : int(f["rc0"]) // 4 - 3] = True
    red_q = (sR - sG > 30.0) & (sR - sB > 35.0) & sint
    red_sq_q = ndimage.binary_opening(red_q, structure=np.ones((8, 8), dtype=bool))
    orange_q = (sR - sG > 35.0) & (sG - sB > 15.0) & (sR > 110.0) & sint
    orange_clean = ndimage.binary_opening(orange_q, structure=np.ones((4, 4), dtype=bool))

    lbl, _ = ndimage.label(red_sq_q)
    sq_centers: list[tuple[int, int]] = []
    for i, sl in enumerate(ndimage.find_objects(lbl), 1):
        m = lbl[sl] == i
        if float((m & orange_clean[sl]).sum()) / float(m.sum()) < 0.25:
            continue
        cy_q, cx_q = ndimage.center_of_mass(m)
        sq_centers.append((int(round((cx_q + sl[1].start) * 4.0 + 2.0)), int(round((cy_q + sl[0].start) * 4.0 + 2.0))))

    red_full = (R - G > 25.0) & (R - B > 25.0)
    red_soft = (R - G > 16.0) & (R - B > 16.0)
    any_ink = gray < 195.0
    hline_r, hline_r0, hline_r1 = find_hline_row_red(image, f)

    r_min_int = int(f["tr1"]) + 10
    r_max_int = int(f["br0"]) - 30
    center_err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0
    cap_err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0

    readings: list[dict] = []
    for c0, r0 in sorted(sq_centers, key=lambda t: t[0]):
        # Refine square center from full-resolution 48x48 red square body
        patch = red_full[r0 - 26 : r0 + 27, c0 - 26 : c0 + 27]
        rw = patch.sum(axis=1) >= 22
        cw = patch.sum(axis=0) >= 22
        body = patch & rw[:, None] & cw[None, :]
        if int(body.sum()) > 200:
            by, bx = ndimage.center_of_mass(body)
            cx_f = float(c0 - 26 + bx + 0.5)
            cy_f = float(r0 - 26 + by + 0.5)
        else:
            cx_f, cy_f = float(c0), float(r0)

        c_int = int(round(cx_f))
        r_int = int(round(cy_f))

        # Full contiguous vertical ink extent along the marker column (black stems/circles may overlay red stem)
        ru = r_int
        while ru > r_min_int and any_ink[ru - 1, c_int - 6 : c_int + 7].any():
            ru -= 1
        rd = r_int
        while rd < r_max_int and any_ink[rd + 1, c_int - 6 : c_int + 7].any():
            rd += 1

        # Red T-caps are ~36-44 px wide (protruding past c +/- 11 on both sides) and terminate
        # inside c +/- 32 (distinguishing them from crossing red curves that continue past c +/- 32).
        # The sole exception is when a bottom red cap terminates directly on the solid horizontal line.
        up_caps: list[int] = []
        dn_caps: list[int] = []
        for r in range(ru, rd + 1):
            if abs(r - r_int) < 32:
                continue
            is_hl = hline_r0 - 6 <= r <= hline_r1 + 6
            lo_ok = red_soft[r, max(0, c_int - 21) : max(0, c_int - 10)].any()
            hi_ok = red_soft[r, c_int + 11 : min(image.shape[1], c_int + 22)].any()
            tot = int(red_soft[r, max(0, c_int - 21) : min(image.shape[1], c_int + 22)].sum())
            if not (lo_ok and hi_ok and tot >= 26):
                continue
            out_l = red_soft[r, max(0, c_int - 42) : max(0, c_int - 32)].any()
            out_r = red_soft[r, c_int + 33 : min(image.shape[1], c_int + 43)].any()
            bounded = (not out_l) and (not out_r)
            hl_term = is_hl and (r > r_int) and (not red_soft[hline_r1 + 12 : min(r_max_int, hline_r1 + 60), c_int - 4 : c_int + 5].any())
            if bounded or hl_term:
                (up_caps if r < r_int else dn_caps).append(r)

        us = _group_runs(np.array(up_caps, dtype=int), gap=4)
        ds = _group_runs(np.array(dn_caps, dtype=int), gap=4)

        if us:
            s0, p0 = us[0]
            # For a top cap merged slightly with stem ink right below, take the top 16-px bar center
            top_cap = float((s0 + min(p0, s0 + 15)) / 2.0 + 0.5)
        else:
            top_cap = float(ru + 4.0)

        if ds:
            s1, p1 = ds[-1]
            bot_cap = float((max(s1, p1 - 15) + p1) / 2.0 + 0.5)
        else:
            bot_cap = float(rd - 4.0)

        x_val = float(x_axis.to_value(cx_f))
        y_val, y_lo, y_hi = value_interval(y_axis, cy_f, center_err_px)
        b_hi = float(y_axis.to_value(top_cap))
        b_lo = float(y_axis.to_value(bot_cap))
        b_hi_out = float(y_axis.to_value(top_cap - cap_err_px))
        b_lo_out = float(y_axis.to_value(bot_cap + cap_err_px))

        readings.append({
            "x": x_val,
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "bar_lo": min(b_lo, b_hi),
            "bar_hi": max(b_lo, b_hi),
            "bar_lo_outer": max(0.0, min(b_lo_out, b_hi_out)),
            "bar_hi_outer": max(b_lo_out, b_hi_out),
            "px": [cx_f, cy_f],
            "px_lo": top_cap,
            "px_hi": bot_cap,
            "note": f"red square (2*kon/D={y_val:.4f}, caps [{b_lo:.4f},{b_hi:.4f}]); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })

    readings.sort(key=lambda r: r["x"])
    return readings
