"""Shared calibrated reader for Fig. 1E (FRAP recovery curves for VASP vs VASP+Arp2/3+VCA)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C

SERIES_RGB = {
    "vasp": np.array([184.0, 112.0, 192.0]),
    "composite": np.array([0.0, 176.0, 80.0]),
}


def read_curve(image: np.ndarray, at: list[float], series: str) -> list[dict]:
    dark, gray, hl, vl = C.find_l_spines(image, thresh=120.0, row_frac=0.45, col_frac=0.45)
    vl0, hb = vl[0], hl[-1]

    yt_cand = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 20), "y", side=-1, band=(5, 35), frac=0.45, max_width=24)
    yt = C.match_ticks_regular(yt_cand, 7, tol=8.0)
    ya = C.verify_axis("linear", list(zip(yt, [120.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0])), max_residual_px=6.0)

    xt_cand = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 20, image.shape[1] - 10), "x", side=+1, band=(5, 35), frac=0.45, max_width=24)
    xt = C.match_ticks_regular(xt_cand, 6, tol=8.0)
    # Calibrate on the first 5 regular ticks (0, 20, 40, 60, 80 s); the 6th equally spaced tick sits at 100 s (misprinted as 120)
    xa = C.verify_axis("linear", list(zip(xt[:5], [0.0, 20.0, 40.0, 60.0, 80.0])), max_residual_px=6.0)

    rgb = image.astype(float)
    mask = np.linalg.norm(rgb - SERIES_RGB[series][None, None, :], axis=2) < 55.0
    # Exclude legend entries in bottom-right interior (t > 35 s, y < 70%)
    r_leg = int(ya.to_pixel(70.0))
    c_leg = int(xa.to_pixel(35.0))
    mask[r_leg:, c_leg:] = False

    r_min = int(yt[0])
    r_max = int(hb[2]) + 15
    out = []
    for x in at:
        out.append(
            C.sample_solid_curve_y(
                mask,
                xa,
                ya,
                float(x),
                r_min,
                r_max,
                search_dx=4,
                note_prefix=f"Fig1E {series}, ",
            )
        )
    return out
