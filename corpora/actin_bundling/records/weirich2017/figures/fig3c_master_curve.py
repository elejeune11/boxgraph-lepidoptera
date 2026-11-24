"""Calibrated reader for Fig. 3C main panel master exponential decay curve."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, ht, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=+1, band=(8, 45))
    xt = C.match_ticks_regular(xt, 5, tol=4.0)
    yt = C.find_ticks(dark, vr0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=-1, band=(8, 45))
    yt = C.match_ticks_regular(yt, 6, tol=4.0)

    xa = C.verify_axis("linear", [(vl0[0], 0.0)] + list(zip(xt, [2.0, 4.0, 6.0, 8.0, 10.0])) + [(vr0[0], 12.0)], 1.5)
    ya = C.verify_axis("linear", list(zip(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])) + [(hb[0], -0.2)], 2.0)
    return xa, ya, int(ht[2]) + 55, int(hb[1]) - 55, int(vl0[2]) + 10, int(vr0[1]) - 55


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    # Pure neutral black pixels isolate the solid exponential curve from all colored scatter symbols
    blk = (r < 40) & (g < 40) & (b < 40) & (np.abs(r - g) < 12) & (np.abs(r - b) < 12)
    mask = np.zeros_like(blk)
    mask[r0:r1, c0:c1] = blk[r0:r1, c0:c1]
    # Mask out the entire inset region (x > 2.6, y > 0.14)
    mask[: int(ya.to_pixel(0.14)), int(xa.to_pixel(2.6)) :] = False
    # Exclude left-spine inward ticks (x < 0.16)
    mask[:, : int(xa.to_pixel(0.16))] = False

    return [
        C.sample_solid_curve_y(mask, xa, ya, float(x), r0, r1, search_dx=6, note_prefix="Fig3C master exp: ")
        for x in at
    ]
