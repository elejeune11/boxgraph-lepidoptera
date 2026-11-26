"""Calibrated reader for Fig. 1(e) (p. 2): measured fascin-actin bundle overtwist vs R (5 red open circles)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import black_mask, calibrate_fig1e, value_interval


def read(image, at):
    x_axis, y_axis, (cL, rT, cR, rB) = calibrate_fig1e(image)
    r_ch = image[:, :, 0].astype(float)
    g_ch = image[:, :, 1].astype(float)
    b_ch = image[:, :, 2].astype(float)
    red = (r_ch > 140.0) & (r_ch - g_ch > 60.0) & (r_ch - b_ch > 60.0)
    blk = black_mask(image, gray_max=80.0, sat_max=35.0)

    # Extend bottom margin by +70 px to capture the R = 0.05 circle centered on the y = 0 axis
    H, W = red.shape
    interior = np.zeros_like(red)
    interior[rT + 20 : min(H, rB + 70), cL + 20 : cR - 20] = red[rT + 20 : min(H, rB + 70), cL + 20 : cR - 20]
    # Exclude legend region (R > 0.23, overtwist < 0.65 deg)
    c_leg = int(round(x_axis.to_pixel(0.23)))
    r_leg = int(round(y_axis.to_pixel(0.65)))
    interior[r_leg:, c_leg:] = False

    lab, n_comp = ndimage.label(interior)
    cal_y = max(y_axis.residual_px, 1.0)
    center_err_px = 2.5 + cal_y
    cap_err_px = 2.0 + cal_y

    readings = []
    for i in range(1, n_comp + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 1000:
            continue
        c_c = 0.5 * (int(cc.min()) + int(cc.max())) + 0.5
        r_c = 0.5 * (int(rr.min()) + int(rr.max())) + 0.5
        r_m_top, r_m_bot = int(rr.min()), int(rr.max())
        c_int = int(round(c_c))

        # Trace black vertical error-bar stems above and below the red circle
        has_up = blk[max(rT + 15, r_m_top - 15) : max(rT + 15, r_m_top - 3), c_int - 12 : c_int + 13].mean() > 0.30
        has_dn = (r_m_bot + 15 < rB) and (blk[r_m_bot + 3 : r_m_bot + 15, c_int - 12 : c_int + 13].mean() > 0.30)

        def trace_stem(r_start: int, direction: int) -> tuple[float, str]:
            r_cur = r_start
            hit_block = False
            hit_axis = False
            while rT + 15 < r_cur < rB - 5 and blk[r_cur, c_int - 12 : c_int + 13].any():
                # Detect entering a filled black espin square (width >= 94 px over >= 25 contiguous rows)
                if direction > 0 and blk[r_cur : min(rB, r_cur + 25), c_int - 55 : c_int + 56].sum(axis=1).min() >= 90:
                    hit_block = True
                    break
                r_cur += direction
            if direction > 0 and r_cur >= rB - 6:
                hit_axis = True
            if hit_block:
                return float("nan"), "occluded_block"
            if hit_axis:
                # Cap merges into horizontal x-axis at y = 0
                return float(y_axis.to_pixel(0.0)), "axis_floor"
            r_end = r_cur - direction
            lo_r = min(r_end, r_end - direction * 25)
            hi_r = max(r_end, r_end - direction * 25)
            cap_rs = [r for r in range(max(rT + 10, lo_r), min(rB - 5, hi_r + 1))
                      if 38 <= blk[r, c_int - 55 : c_int + 56].sum() <= 92]
            if cap_rs:
                return float(np.mean(cap_rs) + 0.5), "cap"
            return float(r_end + 0.5), "stem_end"

        x_val = float(x_axis.to_value(c_c))
        y_val, y_lo, y_hi = value_interval(y_axis, r_c, center_err_px)

        if not has_up and not has_dn:
            top_px, bot_px = float(r_m_top + 0.5), float(r_m_bot + 0.5)
            bar_hi = float(y_axis.to_value(top_px))
            bar_lo = float(y_axis.to_value(bot_px))
            bar_hi_outer = float(y_axis.to_value(top_px - cap_err_px))
            bar_lo_outer = float(y_axis.to_value(bot_px + cap_err_px))
            b_note = "no error bar protruding outside marker: bounded by marker vertical extent"
        else:
            top_px, up_kind = trace_stem(r_m_top - 4, -1) if has_up else (float(r_m_top + 0.5), "marker")
            bot_px, dn_kind = trace_stem(r_m_bot + 4, +1) if has_dn else (float(r_m_bot + 0.5), "marker")
            bar_hi = float(y_axis.to_value(top_px))
            bar_hi_outer = float(y_axis.to_value(top_px - cap_err_px))
            if dn_kind == "occluded_block":
                # Lower error-bar cap hidden inside black espin square at R = 0.33; mirror upper half-bar
                bot_px = 2.0 * r_c - top_px
                bar_lo = 2.0 * y_val - bar_hi
                bar_lo_outer = 2.0 * y_val - bar_hi_outer
                b_note = f"upper cap {top_px:.1f}px; lower cap occluded by espin square (mirrored to {bot_px:.1f}px)"
            else:
                bar_lo = float(y_axis.to_value(bot_px))
                bar_lo_outer = float(y_axis.to_value(bot_px + cap_err_px))
                b_note = f"error bar ({up_kind}={top_px:.1f}px, {dn_kind}={bot_px:.1f}px)"

        readings.append({
            "x": x_val,
            "y": y_val,
            "lo": y_lo,
            "hi": y_hi,
            "bar_lo": min(bar_lo, bar_hi),
            "bar_hi": max(bar_lo, bar_hi),
            "bar_lo_outer": min(bar_lo_outer, bar_hi_outer),
            "bar_hi_outer": max(bar_lo_outer, bar_hi_outer),
            "px": [c_c, r_c],
            "px_lo": top_px,
            "px_hi": bot_px,
            "note": f"fascin circle {len(rr)}px; {b_note}; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px",
        })

    readings.sort(key=lambda d: d["x"])
    return readings
