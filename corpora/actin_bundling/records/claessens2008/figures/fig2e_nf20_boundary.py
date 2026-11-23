"""Calibrated reader for Fig. 2E (p. 2): thick dashed phase-boundary line n_f = 20."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2e import calibrate_fig2e
from _calib import value_interval


def read(image, at):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig2e(image)
    rgb = image.astype(float)
    spread = rgb.max(axis=2) - rgb.min(axis=2)
    black = (rgb.mean(axis=2) < 55.0) & (spread < 25.0)

    interior = np.zeros_like(black)
    interior[int(r_top) + 35 : int(r_bot) - 35, int(c_left) + 35 : int(c_right) - 35] = black[
        int(r_top) + 35 : int(r_bot) - 35, int(c_left) + 35 : int(c_right) - 35
    ]

    # Fit straight line col(row) through the pure-black dash pixels (steep slope ~ -2 on log-log)
    # Use component centroids of clean unobstructed dashes (area >= 400 px)
    lab, n = ndimage.label(interior)
    dc, dr = [], []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) >= 400:
            dc.append(float(cc.mean() + 0.5))
            dr.append(float(rr.mean() + 0.5))
    dc_a, dr_a = np.array(dc), np.array(dr)
    s_cr, i_cr = np.polyfit(dr_a, dc_a, 1)
    inliers = np.abs(dc_a - (s_cr * dr_a + i_cr)) < 12.0
    s_cr, i_cr = np.polyfit(dr_a[inliers], dc_a[inliers], 1)
    fit_res_col = float(np.max(np.abs(dc_a[inliers] - (s_cr * dr_a[inliers] + i_cr))))

    # Convert col-uncertainty to vertical row-uncertainty via 1/|s_cr| (dash thickness ~24 px perpendicular)
    slope_rc = 1.0 / s_cr
    # Perpendicular half-width ~12 px -> vertical half-width = 12 * sqrt(1 + slope_rc^2)
    stroke_half_r = 12.0 * np.sqrt(1.0 + slope_rc * slope_rc)
    half_px = stroke_half_r + abs(slope_rc) * (fit_res_col + x_axis.residual_px) + y_axis.residual_px
    exp_eff = float((y_axis.b / x_axis.b) / s_cr)

    readings = []
    for dd in at:
        c = float(x_axis.to_pixel(dd))
        r = float((c - i_cr) / s_cr)
        y, lo, hi = value_interval(y_axis, r, half_px)
        readings.append({
            "x": float(dd),
            "y": y,
            "lo": lo,
            "hi": hi,
            "px": [c, r],
            "px_lo": r - half_px,
            "px_hi": r + half_px,
            "note": f"dashed n_f=20 line ({int(inliers.sum())} dashes, slope={exp_eff:.3f}, c_a*D_d^2={y*dd*dd:.3f}); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
