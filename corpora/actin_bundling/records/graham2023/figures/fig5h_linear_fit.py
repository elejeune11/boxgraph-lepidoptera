"""Calibrated reader for Main Fig. 5h linear bundle elongation dashed fit line."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.4)
    hb = hl[-1]
    vl0 = vl[0]
    xt = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 10, image.shape[1]), "x", side=+1, band=(6, 30))
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 10), "y", side=-1, band=(6, 30))
    xt = C.match_ticks_regular(xt, 5, tol=5.0)
    yt = C.match_ticks_regular(yt, 7, tol=5.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 50.0, 100.0, 150.0, 200.0])), max_residual_px=3.0)
    ya = C.verify_axis("linear", list(zip(yt, [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0])), max_residual_px=3.0)

    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    # Dashed black fit line: very dark and unsaturated, distinct from colorful semi-transparent traces
    mask = (gray < 55.0) & (sat < 35.0)
    # Exclude axes and upper-left legend ('Measured data' / 'Linear fit')
    mask[int(hb[1]) - 6 :, :] = False
    mask[:, : int(vl0[2]) + 8] = False
    mask[: int(ya.to_pixel(3.5)), : int(xa.to_pixel(150.0))] = False

    row_min = int(yt[0])
    row_max = int(hb[1]) - 6
    return [
        C.sample_dashed_curve_y(mask, xa, ya, x, row_min, row_max, max_search_dx=120, note_prefix="Fig5h linear fit: ")
        for x in at
    ]
