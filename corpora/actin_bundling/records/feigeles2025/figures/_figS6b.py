"""Shared calibrated reader for Fig. S6B (free-solution coalescence length vs time: 12 black dots + red curve)."""

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
    # Calibrate x on the top spine's 7 evenly-spaced ticks (0.2..0.8; tick 0.8 sits 9.5 px inside right spine)
    # matching the 0.050 s data sampling grid to < 0.4 px (bottom ticks are slightly stretched in the PDF graphic)
    xt = C.match_ticks_regular(
        C.find_ticks(dark, ht, (int(vl0[2]) + 1, int(vr0[1]) - 1), "x", side=+1, band=(8, 35)), 7
    )
    yt = C.match_ticks_regular(
        C.find_ticks(dark, vr0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=-1, band=(8, 35)), 5
    )
    xa = C.verify_axis(
        "linear", [(vl0[0], 0.1)] + list(zip(xt, [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])), 2.0
    )
    ya = C.verify_axis(
        "linear", [(ht[0], 7.5)] + list(zip(yt, [7.0, 6.5, 6.0, 5.5, 5.0])) + [(hb[0], 4.5)], 2.0
    )
    return xa, ya, ht, hb, vl0, vr0


def read_data(image: np.ndarray) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[int(ht[2]) + 15 : int(hb[1]) - 10, int(vl0[2]) + 15 : int(vr0[1]) - 15] = True
    # Exclude inset droplet cartoon and equation text (x > 0.27 s, y > 5.4 um)
    inb[: int(ya.to_pixel(5.4)), int(xa.to_pixel(0.27)) :] = False

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
                f"FigS6B black circle (area={len(rr)}px), cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
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
        C.sample_curve_local_fit(red, xa, ya, xv, r0, r1, window_px=65, deg=2, note_prefix="FigS6B red exp fit: ")
        for xv in at
    ]
