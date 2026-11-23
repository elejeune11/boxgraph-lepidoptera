"""Shared calibration and marker extraction for Fig. 3D (p. 3): helical overtwist Delta_theta(R*)."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_ticks, value_interval


def calibrate_fig3d(image: np.ndarray):
    h, v = find_frame(image)
    r_top, r_bot = h[0][2], h[-1][2]
    c_left, c_right = v[0][2], v[-1][2]

    xt_all = find_ticks(image, r_bot, int(c_left) + 8, int(c_right) - 8, "x", -1, offset_lo=8, offset_hi=18)
    xt_maj = [p for p, _, d in xt_all if d >= 30]
    if len(xt_maj) != 3:
        raise RuntimeError(f"Fig. 3D expected 3 major x-ticks (0.01, 0.1, 1), found {len(xt_maj)}")
    x_axis = Axis("log10", [(xt_maj[0], 0.01), (xt_maj[1], 0.1), (xt_maj[2], 1.0)])

    yt_l = find_ticks(image, c_left, int(r_top) + 8, int(r_bot) - 8, "y", +1, offset_lo=8, offset_hi=25)
    y_vals = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
    if len(yt_l) != len(y_vals):
        raise RuntimeError(f"Fig. 3D expected {len(y_vals)} y-ticks, found {len(yt_l)}")
    y_axis = Axis("linear", list(zip([p for p, _, _ in yt_l], y_vals)))
    return x_axis, y_axis, (c_left, r_top, c_right, r_bot)


def read_fig3d_series(image: np.ndarray, which: str):
    x_axis, y_axis, (c_left, r_top, c_right, r_bot) = calibrate_fig3d(image)
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]

    roi = np.zeros(image.shape[:2], bool)
    roi[int(r_top) + 15 : int(r_bot) - 15, int(c_left) + 15 : int(c_right) - 15] = True

    if which == "ca0p8":
        # Filled green triangles; plotted centered at geometric centroid (verified to <0.001 deg)
        mask = roi & (G > 110.0) & (R < 70.0) & (B < 70.0)
        lab, n = ndimage.label(mask)
        use_bbox = False
        label_str = "green triangle (c_a=0.8 mg/ml)"
    elif which == "ca3":
        # Filled orange circles
        mask = roi & (R > 200.0) & (G > 110.0) & (B < 70.0)
        lab, n = ndimage.label(mask)
        use_bbox = False
        label_str = "orange circle (c_a=3 mg/ml)"
    elif which == "ca2":
        # Filled red squares; 3 of 4 are partially overlaid by smaller orange circles or green triangle,
        # leaving all 4 outer corners intact -> dilate to join occluded pieces, then use bounding-box center
        mask = roi & (R > 180.0) & (G < 60.0) & (B < 60.0)
        conn = ndimage.binary_dilation(mask, iterations=12)
        lab, n = ndimage.label(conn)
        use_bbox = True
        label_str = "red square (c_a=2 mg/ml)"
    else:
        raise ValueError(f"unknown Fig. 3D series {which!r}")

    half_px = 3.0 + max(y_axis.residual_px, 0.5)
    readings = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero((lab == i) & mask)
        if len(rr) < 200:
            continue
        if use_bbox:
            r_c = float((rr.min() + rr.max()) / 2.0 + 0.5)
            c_c = float((cc.min() + cc.max()) / 2.0 + 0.5)
        else:
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
            "note": f"{label_str}, {len(rr)} px; cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f}) px",
        })
    readings.sort(key=lambda r: r["x"])
    return readings
