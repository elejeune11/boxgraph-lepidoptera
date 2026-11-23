"""Shared calibration and color masks for Fig. 2E (p. 2): droplet confinement phase diagram."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from boxgraph.digitize import Axis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_ticks


def calibrate_fig2e(image: np.ndarray):
    h, v = find_frame(image)
    r_top, r_bot = h[0][2], h[-1][2]
    c_left, c_right = v[0][2], v[-1][2]

    xt_top = find_ticks(image, r_top, int(c_left) + 8, int(c_right) - 8, "x", +1, offset_lo=8, offset_hi=22)
    xt_bot = find_ticks(image, r_bot, int(c_left) + 8, int(c_right) - 8, "x", -1, offset_lo=8, offset_hi=22)
    x_common = [p for p, _, _ in xt_top if any(abs(p - pb) < 2.0 for pb, _, _ in xt_bot)]
    x_vals = [4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    if len(x_common) != len(x_vals):
        raise RuntimeError(f"Fig. 2E expected {len(x_vals)} x-ticks, found {len(x_common)}")
    x_axis = Axis("log10", list(zip(x_common, x_vals)))

    yt_l = find_ticks(image, c_left, int(r_top) + 8, int(r_bot) - 8, "y", +1, offset_lo=8, offset_hi=22)
    y_vals = [0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]
    if len(yt_l) != 9:
        raise RuntimeError(f"Fig. 2E expected 9 interior y-ticks (0.02..0.10), found {len(yt_l)}")
    y_axis = Axis("log10", list(zip([p for p, _, _ in yt_l] + [r_bot], y_vals)))
    return x_axis, y_axis, (c_left, r_top, c_right, r_bot)


def circle_radius_px(image: np.ndarray, y_axis: Axis, c_right: float) -> float:
    """Measure the exact outer radius of open circle markers from isolated dark-blue circles at c_a=0.1."""
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    m_dblue = (B > 100) & (R < 90) & (G < 100)
    r_01 = int(round(y_axis.to_pixel(0.10)))
    c_start = int(round(c_right - 900))
    c_end = int(round(c_right - 80))
    strip = m_dblue[r_01 - 45 : r_01 + 45, c_start:c_end]
    lab, n = ndimage.label(strip)
    radii = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) > 100:
            radii.append((cc.max() - cc.min() + 1) / 2.0)
    return float(np.median(radii)) if radii else 29.5
