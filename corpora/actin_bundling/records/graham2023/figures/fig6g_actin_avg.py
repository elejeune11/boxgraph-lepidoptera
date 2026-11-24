"""Calibrated reader for Main Fig. 6g average normalized total actin intensity over time (solid black curve)."""

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
    xt = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 10, image.shape[1]), "x", side=+1, band=(6, 25))
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[1]) - 10), "y", side=-1, band=(6, 25), max_width=24)
    xt = C.match_ticks_regular(xt, 5, tol=8.0)
    yt = C.match_ticks_regular(yt, 8, tol=6.0) + [float(hb[0])]
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 50.0, 100.0, 150.0, 200.0])), max_residual_px=5.5)
    ya = C.verify_axis("linear", list(zip(yt, [5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0])), max_residual_px=5.5)

    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    sat = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    # Thick solid black average line, distinct from pastel colored individual curves
    mask = (gray < 45.0) & (sat < 30.0)
    mask[int(hb[1]) - 5 :, :] = False
    mask[:, : int(vl0[2]) + 5] = False
    # Mask upper-left legend
    mask[: int(ya.to_pixel(3.8)), : int(xa.to_pixel(130.0))] = False

    row_min = int(yt[0])
    row_max = int(hb[1]) - 3
    return [
        C.sample_solid_curve_y(mask, xa, ya, x, row_min, row_max, search_dx=3, note_prefix="Fig6g avg intensity: ")
        for x in at
    ]
