"""Shared calibrated reader for Fig. S6C (Delta L = L_i - L_f vs tau: 6 black dots + red linear fit)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, thresh=140.0)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    # Top ticks include 0.6 sitting 12 px inside right spine (exact 382.4 px spacing, residual < 0.4 px)
    xt = C.match_ticks_regular(
        C.find_ticks(dark, ht, (int(vl0[2]) + 1, int(vr0[1]) - 1), "x", side=+1, band=(8, 35)), 5
    )
    # Y-axis has 9 evenly-spaced tick marks (8 equal intervals of 0.2 um from 1.2 bottom spine to 2.8 top spine);
    # the printed text labels "2.8, 2.4, 2, 1.6, 1.2" were placed in a single text box with tighter line spacing
    yt = C.match_ticks_regular(
        C.find_ticks(dark, vr0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=-1, band=(8, 35)), 7
    )
    xa = C.verify_axis(
        "linear", [(vl0[0], 0.1)] + list(zip(xt, [0.2, 0.3, 0.4, 0.5, 0.6])), 2.0
    )
    ya = C.verify_axis(
        "linear",
        [(ht[0], 2.8)] + list(zip(yt, [2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4])) + [(hb[0], 1.2)],
        3.0,
    )
    return xa, ya, ht, hb, vl0, vr0


def read_data(image: np.ndarray) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[int(ht[2]) + 15 : int(hb[1]) - 10, int(vl0[2]) + 15 : int(vr0[1]) - 15] = True

    blk = inb & (R < 60) & (G < 60) & (B < 60)
    lab, n = ndimage.label(blk)
    out = []
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(rr) < 250:
            continue
        c_m = float(cc.mean()) + 0.5
        r_m = float(rr.mean()) + 0.5
        half_px = 2.0 + ya.residual_px
        out.append(
            C.make_marker_reading(
                c_m,
                r_m,
                xa,
                ya,
                half_px,
                f"FigS6C black circle (area={len(rr)}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    out.sort(key=lambda d: d["x"])
    return out


def read_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    r0, r1 = int(ht[2]) + 10, int(hb[1]) - 10
    c0, c1 = int(vl0[2]) + 10, int(vr0[1]) - 10
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    red = inb & (R > 140) & (R - G > 50) & (R - B > 50)
    return [
        C.sample_curve_local_fit(red, xa, ya, xv, r0, r1, window_px=65, deg=1, note_prefix="FigS6C red linear fit: ")
        for xv in at
    ]
