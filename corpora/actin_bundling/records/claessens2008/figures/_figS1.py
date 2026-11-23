"""Shared calibration and curve/marker readers for SI Fig. S1 (p. 6): R* vs filament number n."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_ticks, sample_curve_at, value_interval


def calibrate_figS1(image: np.ndarray):
    h, v = find_frame(image)
    r_top, r_bot = h[0][2], h[-1][2]
    c_left, c_right = v[0][2], v[-1][2]

    xt = find_ticks(image, r_bot, int(c_left) + 8, int(c_right) - 8, "x", -1, offset_lo=8, offset_hi=25)
    x_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    if len(xt) != 17:
        raise RuntimeError(f"Fig. S1 expected 17 interior x-ticks (2..90), found {len(xt)}")
    x_axis = Axis("log10", list(zip([c_left] + [p for p, _, _ in xt] + [c_right], x_vals)))

    yt = find_ticks(image, c_left, int(r_top) + 8, int(r_bot) - 8, "y", +1, offset_lo=8, offset_hi=25)
    y_vals = [
        1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
        0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01,
        0.009, 0.008, 0.007, 0.006,
    ]
    if len(yt) != len(y_vals):
        raise RuntimeError(f"Fig. S1 expected {len(y_vals)} y-ticks, found {len(yt)}")
    y_axis = Axis("log10", list(zip([p for p, _, _ in yt], y_vals)))
    return x_axis, y_axis, (c_left, r_top, c_right, r_bot)


def read_figS1_curve(image: np.ndarray, at: list[float], which: str):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_figS1(image)
    r0, r1 = int(r_top) + 10, int(r_bot) - 10
    c0, c1 = int(c_left) + 10, int(c_right) - 10
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    box = np.zeros(image.shape[:2], bool)
    box[r0:r1, c0:c1] = True

    if which == "blue_surface":
        mask = box & (B > 110.0) & (R < 80.0) & (G < 80.0)
        return [
            sample_curve_at(mask, x_axis, y_axis, nv, r0, r1, window=2, extra_px=1.0, note_prefix="blue surface-included; ")
            for nv in at
        ]
    if which == "green_inside":
        mask = box & (G > 100.0) & (R < 80.0) & (B < 80.0)
        return [
            sample_curve_at(mask, x_axis, y_axis, nv, r0, r1, window=2, extra_px=1.0, note_prefix="green inside-only; ")
            for nv in at
        ]
    raise ValueError(f"unknown Fig. S1 curve {which!r}")


def read_figS1_experiment(image: np.ndarray):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_figS1(image)
    r0, r1 = int(r_top) + 10, int(r_bot) - 10
    c0, c1 = int(c_left) + 10, int(c_right) - 10
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    box = np.zeros(image.shape[:2], bool)
    box[r0:r1, c0:c1] = True
    red = box & (R > 180.0) & (G < 70.0) & (B < 70.0)

    # Experimental filled red circles have diameter ~74 px and no surrounding black border (unlike schematic hexagons)
    yy, xx = np.ogrid[-15:16, -15:16]
    blobs = ndimage.binary_opening(red, structure=(xx * xx + yy * yy) <= 15 * 15)
    lab, n_b = ndimage.label(blobs)

    half_px = 3.0 + max(y_axis.residual_px, 0.5)
    readings = []
    for i in range(1, n_b + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 1500:
            continue
        pr0, pr1 = max(0, int(rr.min()) - 12), min(image.shape[0], int(rr.max()) + 12)
        pc0, pc1 = max(0, int(cc.min()) - 12), min(image.shape[1], int(cc.max()) + 12)
        patch = rgb[pr0:pr1, pc0:pc1]
        black_frac = ((patch.mean(axis=2) < 60.0) & ((patch.max(axis=2) - patch.min(axis=2)) < 40.0)).mean()
        if black_frac >= 0.05:
            continue  # schematic hexagon circle with black outline

        r_c = float(rr.mean() + 0.5)
        c_c = float(cc.mean() + 0.5)
        x_val = float(x_axis.to_value(c_c))
        y_val, lo, hi = value_interval(y_axis, r_c, half_px)
        readings.append({
            "x": x_val,
            "y": y_val,
            "lo": lo,
            "hi": hi,
            "px": [c_c, r_c],
            "px_lo": r_c - half_px,
            "px_hi": r_c + half_px,
            "note": f"red experimental circle ({len(rr)} px); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    readings.sort(key=lambda r: r["y"])
    return readings


def read_figS1_powerlaw(image: np.ndarray, at: list[float]):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_figS1(image)
    r0, r1 = int(r_top) + 10, int(r_bot) - 10
    c0, c1 = int(c_left) + 10, int(c_right) - 10
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    box = np.zeros(image.shape[:2], bool)
    box[r0:r1, c0:c1] = True
    red = box & (R > 180.0) & (G < 70.0) & (B < 70.0)

    # Remove all red filled circles (experimental and schematic) to isolate the thin red dashed line
    yy, xx = np.ogrid[-12:13, -12:13]
    red_thick = ndimage.binary_dilation(ndimage.binary_opening(red, structure=(xx * xx + yy * yy) <= 12 * 12), iterations=16)
    red_dash = red & ~red_thick

    # Restrict to dashed-line bounding region n in [1.5, 20], R* in [0.007, 0.36]
    cc0, cc1 = int(x_axis.to_pixel(1.5)), int(x_axis.to_pixel(20.0))
    rr0, rr1 = int(y_axis.to_pixel(0.36)), int(y_axis.to_pixel(0.007))
    roi = np.zeros_like(red_dash)
    roi[rr0:rr1, cc0:cc1] = red_dash[rr0:rr1, cc0:cc1]

    lab, n_d = ndimage.label(roi)
    dc, dr = [], []
    for i in range(1, n_d + 1):
        rr, cc = np.nonzero(lab == i)
        if 80 <= len(rr) <= 2000:
            dc.append(float(cc.mean() + 0.5))
            dr.append(float(rr.mean() + 0.5))
    dc_a, dr_a = np.array(dc), np.array(dr)
    slope, intercept = np.polyfit(dc_a, dr_a, 1)
    inliers = np.abs(dr_a - (slope * dc_a + intercept)) < 12.0
    slope, intercept = np.polyfit(dc_a[inliers], dr_a[inliers], 1)
    fit_max_res = float(np.max(np.abs(dr_a[inliers] - (slope * dc_a[inliers] + intercept))))

    half_px = 6.0 + fit_max_res + max(y_axis.residual_px, 0.5) + abs(slope) * max(x_axis.residual_px, 0.5)
    exp_eff = float((y_axis.b / x_axis.b) * slope)

    readings = []
    for nv in at:
        c = float(x_axis.to_pixel(nv))
        r = float(slope * c + intercept)
        y, lo, hi = value_interval(y_axis, r, half_px)
        readings.append({
            "x": float(nv),
            "y": y,
            "lo": lo,
            "hi": hi,
            "px": [c, r],
            "px_lo": r - half_px,
            "px_hi": r + half_px,
            "note": f"red dashed power-law line ({int(inliers.sum())} dashes, slope={exp_eff:.4f}, fit_max_res={fit_max_res:.2f} px); cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    return readings
