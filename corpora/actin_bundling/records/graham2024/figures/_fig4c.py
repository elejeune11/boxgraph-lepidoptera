"""Shared calibrated reader for Fig. 4C bottom 1D actin intensity line profiles."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read_profile(image: np.ndarray, at: list[float], label: str) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=110.0, row_frac=0.4, col_frac=0.35)
    vl0 = vl[0]
    y_end = np.where(dark[:, int(round(vl0[0]))])[0].max()
    hb = min(hl, key=lambda h: abs(h[0] - y_end))

    yt_c = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 25), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_c, 4, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [1.2, 0.8, 0.4, 0.0])), max_residual_px=4.0)

    xt_c = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 25, image.shape[1] - 5), "x", side=+1, band=(5, 35), frac=0.45, max_width=24)
    xt = C.match_ticks_regular(xt_c, 5, tol=10.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 1.0, 2.0, 3.0, 4.0])), max_residual_px=6.0)

    rgb = image.astype(float)
    grn = (rgb[:, :, 1] > rgb[:, :, 0] + 35.0) & (rgb[:, :, 1] > rgb[:, :, 2] + 20.0)
    r_min = int(yt[0]) - 10
    r_max = int(hb[2]) + 15
    return [
        C.sample_solid_curve_y(grn, xa, ya, float(x), r_min, r_max, search_dx=6, note_prefix=f"Fig4C {label}, ")
        for x in at
    ]
