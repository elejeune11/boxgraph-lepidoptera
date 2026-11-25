"""Shared calibrated reader for Fig. 1E (FUS blue and actin green transverse intensity profiles)."""

from __future__ import annotations

import numpy as np
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, thresh=140.0)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt_t = C.match_ticks_regular(
        C.find_ticks(dark, ht, (int(vl0[2]) + 15, int(vr0[1]) - 15), "x", side=+1, band=(4, 18)), 4
    )
    yt_l = C.match_ticks_regular(
        C.find_ticks(dark, vl0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=+1, band=(4, 18)), 5
    )
    xa = C.verify_axis("linear", [(vl0[0], 0.0)] + list(zip(xt_t, [2.0, 4.0, 6.0, 8.0])) + [(vr0[0], 10.0)], 2.5)
    ya = C.verify_axis(
        "linear", [(ht[0], 800.0)] + list(zip(yt_l, [700.0, 600.0, 500.0, 400.0, 300.0])) + [(hb[0], 200.0)], 2.5
    )
    return xa, ya, int(ht[2]) + 8, int(hb[1]) - 8, int(vl0[2]) + 8, int(vr0[1]) - 8


def read_curve(image: np.ndarray, at: list[float], channel: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True

    if channel == "fus":
        # Bright cyan-blue FUS curve
        mask = inb & (B - R > 60) & (B > 140)
        prefix = "Fig1E FUS (cyan): "
    elif channel == "actin":
        # Forest-green actin curve (exclude cyan pixels where B > G - 10)
        mask = inb & (G - R > 30) & (G - B > 15) & (G > 70)
        prefix = "Fig1E actin (green): "
    else:
        raise ValueError(channel)

    return [
        C.sample_solid_curve_y(mask, xa, ya, xv, r0, r1, search_dx=10, note_prefix=prefix)
        for xv in at
    ]
