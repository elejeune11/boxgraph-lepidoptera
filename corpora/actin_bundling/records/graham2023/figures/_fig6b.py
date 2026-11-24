"""Shared calibrated reader for Main Fig. 6b FRAP average recovery curves (3% PEG red, 10% PEG blue)."""

from __future__ import annotations

import numpy as np
import _calib as C


def read_frap(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.3, col_frac=0.4)
    vl0 = vl[0]
    yt_all = C.find_outward_ticks(dark, vl0, (0, image.shape[0]), "y", side=-1, band=(6, 25))
    # Top 6 ticks (120..20) are completely unobstructed by the cyan curve on y=0 (exact 0.00 px linear fit)
    yt = C.match_ticks_regular(yt_all, 6, tol=4.0)
    ya = C.verify_axis("linear", list(zip(yt, [120.0, 100.0, 80.0, 60.0, 40.0, 20.0])), max_residual_px=2.0)
    r_zero = int(round(ya.to_pixel(0.0)))
    xt = C.find_outward_ticks(dark, (float(r_zero), r_zero - 2, r_zero + 2), (int(vl0[1]) - 10, image.shape[1]), "x", side=+1, band=(10, 35))
    xt = C.match_ticks_regular(xt, 4, tol=5.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 100.0, 200.0, 300.0])), max_residual_px=3.0)

    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    if which == "3peg":
        # Bright red curve, distinct from grey individual traces
        mask = (r > 170.0) & (r - g > 90.0) & (r - b > 90.0)
        label = "Fig6b 3% PEG: "
    elif which == "10peg":
        # Bright cyan-blue curve
        mask = (b > 170.0) & (b - r > 70.0) & (g > 130.0)
        label = "Fig6b 10% PEG: "
    else:
        raise ValueError(which)

    # Mask upper legend and micrograph inset region (y > 70% for x > 20 s)
    mask[: int(ya.to_pixel(65.0)), int(xa.to_pixel(25.0)) :] = False

    row_min = int(yt[0]) - 20
    row_max = r_zero + 25
    return [
        C.sample_solid_curve_y(mask, xa, ya, x, row_min, row_max, search_dx=3, note_prefix=label)
        for x in at
    ]
