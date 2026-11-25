"""Calibrated reader for Fig. S7 (histogram of capillary pressure for 28 FUS condensates, with SD error bars)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_box_spines(image, thresh=140.0)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.match_ticks_regular(
        C.find_ticks(dark, ht, (int(vl0[2]) + 15, int(vr0[1]) - 15), "x", side=+1, band=(8, 35)), 5
    )
    yt = C.match_ticks_regular(
        C.find_ticks(dark, vr0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=-1, band=(8, 35)), 5
    )
    xa = C.verify_axis(
        "linear", [(vl0[0], 10.0)] + list(zip(xt, [20.0, 30.0, 40.0, 50.0, 60.0])) + [(vr0[0], 70.0)], 3.0
    )
    ya = C.verify_axis(
        "linear", [(ht[0], 12.0)] + list(zip(yt, [10.0, 8.0, 6.0, 4.0, 2.0])) + [(hb[0], 0.0)], 3.0
    )
    return C.read_hist_bars_red_err(image, xa, ya, ht, hb, vl0, vr0, note_prefix="FigS7 hist: ")
