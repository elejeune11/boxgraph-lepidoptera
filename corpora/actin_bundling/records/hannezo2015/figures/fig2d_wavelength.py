"""Calibrated reader for Main Fig. 2D right (p. 3): dimensionless wavelength lambda_c/l_h vs turnover phi."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import find_frame, find_inward_ticks, match_ticks_regular, verify_axis, sample_solid_curve_y


def read(image, at):
    dark, _, hl, vl = find_frame(image, row_frac=0.55, col_frac=0.55, thresh=140.0)
    h_top, h_bot = hl[0], hl[-1]
    v_left, v_right = vl[0], vl[-1]
    xt_top, _, yt_left, _ = find_inward_ticks(dark, h_top, h_bot, v_left, v_right, band=(4, 15), margin=12, max_width=10)

    # Note: the 3 interior ticks at 4, 6, 8 in this vector figure are shifted ~26.5 px left of their
    # printed labels '4', '6', '8' (an Illustrator grouping artifact), whereas left border (0),
    # first interior tick (2), and right border (10) are collinear to 0.15 px and align with all labels.
    x2 = min(xt_top)
    x_axis = verify_axis("linear", [(v_left[0], 0.0), (x2, 2.0), (v_right[0], 10.0)], max_residual_px=1.5)

    # Y ticks: top border = 15, 7 interior ticks = 14, 12, 10, 8, 6, 4, 2, bottom border = 0
    yt_int = match_ticks_regular(yt_left, 7, tol=2.5)
    y_axis = verify_axis(
        "linear",
        [(h_top[0], 15.0)] + list(zip(yt_int, [14.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0])) + [(h_bot[0], 0.0)],
        max_residual_px=1.5,
    )

    rgb = image.astype(float)
    red = (rgb[:, :, 0] > 180) & (rgb[:, :, 1] < 60) & (rgb[:, :, 2] < 60)
    # Exclude legend line in top right
    red[: int(h_top[2] + 150), int(x_axis.to_pixel(6.5)) :] = False

    r_min, r_max = int(h_top[2]) + 10, int(h_bot[1]) - 5
    return [sample_solid_curve_y(red, x_axis, y_axis, x, r_min, r_max) for x in at]
