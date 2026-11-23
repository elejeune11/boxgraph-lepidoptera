"""Calibrated reader for Fig. 1C (p. 2): dashed power-law line D ~ (R*)^0.3."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig1c_diameter_data import calibrate_fig1c
from _calib import value_interval


def read(image, at):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig1c(image)
    gray = image.astype(float).mean(axis=2)
    dark = gray < 100.0

    interior = np.zeros_like(dark)
    r0, r1 = int(r_top) + 48, int(r_bot) - 48
    c0, c1 = int(c_left) + 48, int(c_right) - 15
    interior[r0:r1, c0:c1] = dark[r0:r1, c0:c1]
    r_ins = int(round(y_axis.to_pixel(15.5)))
    c_ins = int(round(x_axis.to_pixel(0.095)))
    interior[r_ins:, c_ins:] = False

    # Mask out the 8 circular markers and their error-bar columns
    yy, xx = np.ogrid[-18:19, -18:19]
    markers = ndimage.binary_opening(interior, structure=(xx * xx + yy * yy) <= 18 * 18)
    lab_m, nm = ndimage.label(markers)
    clean_dash = interior.copy()
    for i in range(1, nm + 1):
        rr, cc = np.nonzero(lab_m == i)
        if len(rr) >= 1000:
            c_c = int(round(cc.mean()))
            clean_dash[:, max(0, c_c - 48) : min(clean_dash.shape[1], c_c + 48)] = False

    # Fit straight line through isolated dash-segment centroids (power law is straight on log-log axes)
    lab_d, nd = ndimage.label(clean_dash)
    dash_c, dash_r = [], []
    for i in range(1, nd + 1):
        rr, cc = np.nonzero(lab_d == i)
        if 150 <= len(rr) <= 2000:
            dash_c.append(float(cc.mean() + 0.5))
            dash_r.append(float(rr.mean() + 0.5))
    dc, dr = np.array(dash_c), np.array(dash_r)
    slope, intercept = np.polyfit(dc, dr, 1)
    res_dash = np.abs(dr - (slope * dc + intercept))
    inliers = res_dash < 5.0
    slope, intercept = np.polyfit(dc[inliers], dr[inliers], 1)
    max_fit_res = float(np.max(np.abs(dr[inliers] - (slope * dc[inliers] + intercept))))

    # Vertical stroke half-width of dashes (~5 px) + line fit residual + calibration residuals
    half_px = 5.0 + max_fit_res + max(y_axis.residual_px, 0.5) + abs(slope) * max(x_axis.residual_px, 0.5)
    exp_eff = float((y_axis.b / x_axis.b) * slope)

    readings = []
    for x_val in at:
        c = float(x_axis.to_pixel(x_val))
        r = float(slope * c + intercept)
        y, lo, hi = value_interval(y_axis, r, half_px)
        readings.append({
            "x": float(x_val),
            "y": y,
            "lo": lo,
            "hi": hi,
            "px": [c, r],
            "px_lo": r - half_px,
            "px_hi": r + half_px,
            "note": f"dashed line ({int(inliers.sum())} dashes, slope={exp_eff:.4f}, fit_max_res={max_fit_res:.2f} px); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
