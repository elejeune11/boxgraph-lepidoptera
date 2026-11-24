"""Calibrated reader for Main Fig. 1f FRAP average recovery curve (red line, n = 9 droplets)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.4)
    hb = hl[-1]
    vl0 = vl[0]  # leftmost vertical spine is the y-axis (near -10 s; vertical red drop is slightly right)
    xt = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 10, image.shape[1]), "x", side=+1, band=(8, 35))
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 10), "y", side=-1, band=(8, 35))
    xt = C.match_ticks_regular(xt, 7, tol=8.0)
    yt = C.match_ticks_regular(yt, 7, tol=8.0)
    xa = C.verify_axis("linear", list(zip(xt, [-10.0, 40.0, 90.0, 140.0, 190.0, 240.0, 290.0])), max_residual_px=6.0)
    ya = C.verify_axis("linear", list(zip(yt, [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0])), max_residual_px=6.0)

    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    # Bright red average curve, completely separable from grey individual trials
    red = (r > 170.0) & (r - g > 90.0) & (r - b > 90.0)

    row_min = int(yt[0]) - 20
    row_max = int(hb[2]) + 30
    return [
        C.sample_solid_curve_y(red, xa, ya, x, row_min, row_max, search_dx=3, note_prefix="Fig1f FRAP avg: ")
        for x in at
    ]
