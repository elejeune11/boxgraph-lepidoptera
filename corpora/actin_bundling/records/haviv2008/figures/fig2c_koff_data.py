"""Calibrated reader for Fig. 2c (p. 4 / 450): fitted 10*k_off/D (divided by 10 -> k_off/D), 5 black circles + ±20% bars."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import NATIVE_QUANT_PX, _group_runs, calibrate_fig2c, value_interval


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    x_axis, y_axis, f = calibrate_fig2c(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    gray = (R + G + B) / 3.0

    # Detect the 5 black hollow circles on 1/4-res via hole-fill + disk opening
    sub = image[::4, ::4].astype(float)
    sR, sG, sB = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    sgray = (sR + sG + sB) / 3.0
    sint = np.zeros(sgray.shape, dtype=bool)
    sint[int(f["tr1"]) // 4 + 3 : int(f["br0"]) // 4 - 3, int(f["lc1"]) // 4 + 3 : int(f["rc0"]) // 4 - 3] = True
    bk_q = (sgray < 160.0) & (sR - sG < 25.0) & sint
    bk_filled = ndimage.binary_fill_holes(bk_q)
    yy, xx = np.ogrid[-4:5, -4:5]
    disk9 = (xx * xx + yy * yy) <= 17
    bk_circ = ndimage.binary_opening(bk_filled, structure=disk9)
    lbl, _ = ndimage.label(bk_circ)

    # Full-resolution dark neutral mask for center refinement & black stem/cap tracing
    bk_neutral = (gray < 155.0) & (R - G < 28.0)
    # Dark cap mask allows slight warm edge tint where a black cap touches an adjacent red square
    bk_cap = (gray < 145.0) & (R - G < 48.0)
    any_ink = gray < 195.0

    r_min_int = int(f["tr1"]) + 10
    r_max_int = int(f["br0"]) - 30  # stay above bottom inward ticks

    center_err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0
    cap_err_px = NATIVE_QUANT_PX + max(y_axis.residual_px, 1.0) + 2.0

    readings: list[dict] = []
    for i, sl in enumerate(ndimage.find_objects(lbl), 1):
        m = lbl[sl] == i
        cy_q, cx_q = ndimage.center_of_mass(m)
        c0 = int(round((cx_q + sl[1].start) * 4.0 + 2.0))
        r0 = int(round((cy_q + sl[0].start) * 4.0 + 2.0))

        # Refine circle center from full-resolution circle ring bounding box
        patch = bk_neutral[r0 - 28 : r0 + 29, c0 - 28 : c0 + 29]
        p_fill = ndimage.binary_fill_holes(patch)
        # Keep rows/cols that belong to the ~40x40 circular body (width >= 20 px)
        rw = p_fill.sum(axis=1) >= 20
        cw = p_fill.sum(axis=0) >= 20
        body = p_fill & rw[:, None] & cw[None, :]
        if int(body.sum()) > 100:
            by, bx = ndimage.center_of_mass(body)
            cx_f = float(c0 - 28 + bx + 0.5)
            cy_f = float(r0 - 28 + by + 0.5)
        else:
            cx_f, cy_f = float(c0), float(r0)

        c_int = int(round(cx_f))
        r_int = int(round(cy_f))

        # Find contiguous vertical stem extent containing r_int (using any_ink so overlapping red markers don't cut stem)
        r_up_lim = r_int
        while r_up_lim > r_min_int and any_ink[r_up_lim - 1, c_int - 6 : c_int + 7].any():
            r_up_lim -= 1
        r_dn_lim = r_int
        while r_dn_lim < r_max_int and any_ink[r_dn_lim + 1, c_int - 6 : c_int + 7].any():
            r_dn_lim += 1

        # Scan candidate horizontal black T-caps along [r_up_lim, r_int - 28] and [r_int + 28, r_dn_lim]
        def _find_black_cap(r_start: int, r_end: int, direction: str) -> float:
            # Find runs where bk_cap extends symmetrically on BOTH left (c-15..c-5) and right (c+6..c+16)
            rows_in = []
            for r in range(max(r_min_int, r_start), min(r_max_int, r_end) + 1):
                left_ok = bk_cap[r, c_int - 15 : c_int - 4].any()
                right_ok = bk_cap[r, c_int + 5 : c_int + 16].any()
                total = int(bk_cap[r, c_int - 15 : c_int + 16].sum())
                # Also exclude diagonal crossing curves by requiring stem is centered and cap is narrow vertically
                if left_ok and right_ok and total >= 18:
                    rows_in.append(r)
            segs = _group_runs(np.array(rows_in, dtype=int), gap=2)
            # Filter out segments that are too tall (marker bodies > 28 px tall) or wide curves far from cap
            cap_segs = [(s, p) for s, p in segs if (p - s + 1) <= 24]
            # For top cap: take the highest black cap whose rows have neutral-black presence (not pure red)
            valid = []
            for s, p in cap_segs:
                # Confirm genuine dark cap (median darkness < 135)
                if float(gray[s : p + 1, c_int - 12 : c_int + 13].min()) < 120.0:
                    # Check there is no dark neutral stem continuing well past the cap (within 25..45 px)
                    if direction == "up":
                        above = bk_neutral[max(r_min_int, s - 40) : max(r_min_int, s - 8), c_int - 4 : c_int + 5]
                        if above.size == 0 or above.mean() < 0.25:
                            valid.append((s, p))
                    else:
                        below = bk_neutral[min(r_max_int, p + 9) : min(r_max_int, p + 41), c_int - 4 : c_int + 5]
                        if below.size == 0 or below.mean() < 0.25:
                            valid.append((s, p))
            if valid:
                chosen = valid[0] if direction == "up" else valid[-1]
                return float((chosen[0] + chosen[1]) / 2.0 + 0.5)
            # Fallback to outermost dark neutral row on the stem
            if direction == "up":
                nz = np.where(bk_cap[r_start : r_end + 1, c_int - 6 : c_int + 7].any(axis=1))[0]
                return float(r_start + nz[0] + 4.0) if len(nz) else float(r_start + 4.0)
            nz = np.where(bk_cap[r_start : r_end + 1, c_int - 6 : c_int + 7].any(axis=1))[0]
            return float(r_start + nz[-1] - 4.0) if len(nz) else float(r_end - 4.0)

        top_cap = _find_black_cap(r_up_lim, r_int - 28, "up")
        bot_cap = _find_black_cap(r_int + 28, r_dn_lim, "down")

        x_val = float(x_axis.to_value(cx_f))
        # Axis displays 10 * k_off / D; divide by 10 to report k_off / D
        y10, y10_lo, y10_hi = value_interval(y_axis, cy_f, center_err_px)
        b_hi10 = float(y_axis.to_value(top_cap))
        b_lo10 = float(y_axis.to_value(bot_cap))
        b_hi10_out = float(y_axis.to_value(top_cap - cap_err_px))
        b_lo10_out = float(y_axis.to_value(bot_cap + cap_err_px))

        readings.append({
            "x": x_val,
            "y": y10 / 10.0,
            "lo": y10_lo / 10.0,
            "hi": y10_hi / 10.0,
            "bar_lo": min(b_lo10, b_hi10) / 10.0,
            "bar_hi": max(b_lo10, b_hi10) / 10.0,
            "bar_lo_outer": max(0.0, min(b_lo10_out, b_hi10_out) / 10.0),
            "bar_hi_outer": max(b_lo10_out, b_hi10_out) / 10.0,
            "px": [cx_f, cy_f],
            "px_lo": top_cap,
            "px_hi": bot_cap,
            "note": f"black circle (10*koff/D={y10:.4f}, caps [{b_lo10:.4f},{b_hi10:.4f}], div by 10); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })

    readings.sort(key=lambda r: r["x"])
    return readings
