"""Calibrated reader for SI Fig. S5 (p. 12): linear stability boundary chi*rho_0 vs xi (log10 x-axis)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import verify_axis, sample_solid_curve_y


def read(image, at):
    rgb = image.astype(float)
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]

    # Detect exact colored phase-diagram grid box (xi in [10^-5, 10^0], chi*rho_0 in [0, 10])
    colored = (np.abs(R - G) > 25.0) | (np.abs(R - B) > 25.0) | (np.abs(G - B) > 25.0)
    # Exclude any rightmost legend sliver if present (grid ends at black outer vertical border)
    H, W = colored.shape
    c_end = W - 1
    # Find rightmost contiguous colored block belonging to the main grid
    col_has = np.where(colored.mean(axis=0) > 0.35)[0]
    c_left = float(col_has.min()) - 0.5
    c_right = float(col_has.max()) + 0.5
    row_has = np.where(colored[:, int(c_left) : int(c_right) + 1].mean(axis=1) > 0.35)[0]
    r_top = float(row_has.min()) - 0.5
    r_bot = float(row_has.max()) + 0.5

    # Verify outer black borders sit immediately around colored box
    gray = rgb.mean(axis=2)
    # Refine box edges to center of outer black border lines
    cl_black = np.where(gray[int(r_top) + 20 : int(r_bot) - 20, max(0, int(c_left) - 10) : int(c_left) + 5].mean(axis=0) < 90)[0] + max(0, int(c_left) - 10)
    cr_black = np.where(gray[int(r_top) + 20 : int(r_bot) - 20, int(c_right) - 4 : min(W, int(c_right) + 10)].mean(axis=0) < 90)[0] + int(c_right) - 4
    rt_black = np.where(gray[max(0, int(r_top) - 10) : int(r_top) + 5, int(c_left) + 20 : int(c_right) - 20].mean(axis=1) < 90)[0] + max(0, int(r_top) - 10)
    rb_black = np.where(gray[int(r_bot) - 4 : min(H, int(r_bot) + 10), int(c_left) + 20 : int(c_right) - 20].mean(axis=1) < 90)[0] + int(r_bot) - 4

    c0 = float(cl_black.mean()) + 0.5 if len(cl_black) else c_left
    c1 = float(cr_black.mean()) + 0.5 if len(cr_black) else c_right
    r0 = float(rt_black.mean()) + 0.5 if len(rt_black) else r_top
    r1 = float(rb_black.mean()) + 0.5 if len(rb_black) else r_bot

    x_axis = verify_axis("log10", [(c0, 1e-5), (c1, 1.0)], max_residual_px=1.5)
    y_axis = verify_axis("linear", [(r0, 10.0), (r1, 0.0)], max_residual_px=1.5)

    # Bright red linear stability curve (distinct from dark-red chaos cells where R ~ 128)
    bright_red = (R > 210.0) & (G < 55.0) & (B < 55.0)
    r_min, r_max = int(r0) + 5, int(r1) - 5
    return [sample_solid_curve_y(bright_red, x_axis, y_axis, x, r_min, r_max, search_dx=3) for x in at]
