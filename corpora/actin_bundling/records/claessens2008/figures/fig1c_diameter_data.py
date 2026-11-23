"""Calibrated reader for Fig. 1C (p. 2): measured bundle diameter D(R*) with error bars."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_ticks, value_interval


def calibrate_fig1c(image: np.ndarray):
    h, v = find_frame(image)
    r_top, r_bot = h[0][2], h[-1][2]
    c_left, c_right = v[0][2], v[-1][2]

    xt_top = find_ticks(image, r_top, int(c_left) + 8, int(c_right) - 8, "x", +1, offset_lo=8, offset_hi=18)
    x_maj = [p for p, _, d in xt_top if d > 30]
    if len(x_maj) != 2:
        raise RuntimeError(f"Fig. 1C expected 2 interior major x-ticks (0.01, 0.1), found {len(x_maj)}")
    x_axis = Axis("log10", [(x_maj[0], 0.01), (x_maj[1], 0.1), (c_right, 1.0)])

    yt_l = find_ticks(image, c_left, int(r_top) + 8, int(r_bot) - 8, "y", +1, offset_lo=8, offset_hi=18)
    yt_r = find_ticks(image, c_right, int(r_top) + 8, int(r_bot) - 8, "y", -1, offset_lo=8, offset_hi=18)
    y_matched = [p for p, _, _ in yt_l if any(abs(p - pr) < 3.0 for pr, _, _ in yt_r)]
    if len(y_matched) != 6:
        raise RuntimeError(f"Fig. 1C expected 6 matched y-ticks (30,20,10,9,8,7), found {len(y_matched)}")
    y_axis = Axis("log10", list(zip(y_matched, [30.0, 20.0, 10.0, 9.0, 8.0, 7.0])))
    return x_axis, y_axis, (c_left, r_top, c_right, r_bot)


def read(image, at):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig1c(image)
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0

    # Restrict to plot interior excluding ticks and the bottom-right TEM micrograph inset
    interior = np.zeros_like(dark)
    r0, r1 = int(r_top) + 48, int(r_bot) - 48
    c0, c1 = int(c_left) + 48, int(c_right) - 15
    interior[r0:r1, c0:c1] = dark[r0:r1, c0:c1]
    # Exclude TEM inset (located in bottom-right: R* > 0.095, D < 15.5 nm)
    r_ins = int(round(y_axis.to_pixel(15.5)))
    c_ins = int(round(x_axis.to_pixel(0.095)))
    interior[r_ins:, c_ins:] = False

    # Separate filled circular markers (diameter ~58 px) from error bars and dashes
    yy, xx = np.ogrid[-18:19, -18:19]
    disk = (xx * xx + yy * yy) <= 18 * 18
    markers = ndimage.binary_opening(interior, structure=disk)
    lab, n_m = ndimage.label(markers)

    cal_y = max(y_axis.residual_px, 0.5)
    cap_err_px = 1.5 + cal_y
    center_err_px = 2.0 + cal_y

    readings = []
    for i in range(1, n_m + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 1000:
            continue
        r_c = float(rr.mean() + 0.5)
        c_c = float(cc.mean() + 0.5)
        c_int = int(round(c_c))

        # Trace vertical error bar stem above and below marker
        r_up = int(rr.min())
        while r_up > int(r_top) + 12 and dark[r_up - 1, c_int - 3 : c_int + 4].any():
            r_up -= 1
        r_dn = int(rr.max())
        while r_dn < int(r_bot) - 12 and dark[r_dn + 1, c_int - 3 : c_int + 4].any():
            r_dn += 1

        # Locate horizontal error-bar caps (width >= 12 px, vs stem width ~5 px)
        cap_up = [r for r in range(r_up, min(r_up + 12, int(rr.min()))) if dark[r, c_int - 15 : c_int + 16].sum() >= 12]
        cap_dn = [r for r in range(max(r_dn - 12, int(rr.max()) + 1), r_dn + 1) if dark[r, c_int - 15 : c_int + 16].sum() >= 12]
        top_cap = float(np.mean(cap_up) + 0.5) if cap_up else float(r_up + 0.5)
        bot_cap = float(np.mean(cap_dn) + 0.5) if cap_dn else float(r_dn + 0.5)

        x_val = float(x_axis.to_value(c_c))
        y_val, y_lo, y_hi = value_interval(y_axis, r_c, center_err_px)
        bar_hi = float(y_axis.to_value(top_cap))
        bar_lo = float(y_axis.to_value(bot_cap))
        bar_hi_outer = float(y_axis.to_value(top_cap - cap_err_px))
        bar_lo_outer = float(y_axis.to_value(bot_cap + cap_err_px))

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
            "px_lo": top_cap,
            "px_hi": bot_cap,
            "note": f"marker {len(rr)} px, error bar caps ({top_cap:.1f}, {bot_cap:.1f}) px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })

    readings.sort(key=lambda r: r["x"])
    return readings
