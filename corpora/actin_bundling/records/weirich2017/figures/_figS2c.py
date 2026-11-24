"""Shared calibrated reader for SI Fig. S2C (rescaled free energy vs alpha)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt_major = C.find_ticks(dark, hb, (int(vl0[2]) + 10, int(vr0[1]) - 10), "x", side=-1, band=(8, 45))
    xt_major = C.match_ticks_regular(xt_major, 6, tol=4.0)
    yt_major = C.find_ticks(dark, vl0, (int(ht[2]) + 10, int(hb[1]) - 10), "y", side=+1, band=(8, 45))
    yt_major = C.match_ticks_regular(yt_major, 3, tol=4.0)

    x_marks = [(vl0[0], 0.0)] + list(zip(xt_major, [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]))
    y_marks = [(ht[0], 6.0)] + list(zip(yt_major, [5.0, 4.0, 3.0])) + [(hb[0], 2.0)]
    xa = C.verify_axis("linear", x_marks, max_residual_px=1.0)
    ya = C.verify_axis("linear", y_marks, max_residual_px=1.0)
    return xa, ya, int(ht[2]) + 4, int(hb[1]) - 4, int(vl0[2]) + 4, int(vr0[1]) - 4


def _curve_mask(image: np.ndarray, which: str, r0: int, r1: int, c0: int, c1: int) -> np.ndarray:
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)
    if which == "w2":
        # Green curve RGB ~(26, 167, 53)
        raw = np.sqrt((r - 26.0) ** 2 + (g - 167.0) ** 2 + (b - 53.0) ** 2) < 70.0
    elif which == "w0p5":
        # Orange/gold curve RGB ~(253, 195, 0)
        raw = np.sqrt((r - 253.0) ** 2 + (g - 195.0) ** 2 + (b - 0.0) ** 2) < 85.0
    elif which == "w0":
        # Blue curve RGB ~(95, 129, 183)
        raw = np.sqrt((r - 95.0) ** 2 + (g - 129.0) ** 2 + (b - 183.0) ** 2) < 70.0
    else:
        raise ValueError(which)

    mask = np.zeros_like(raw)
    mask[r0:r1, c0:c1] = raw[r0:r1, c0:c1]

    # Remove decorative solid tactoid icons (thick filled blobs) and text labels:
    # The thin curve is the single component with the largest horizontal span (> 1500 px).
    lab, n = ndimage.label(mask)
    keep = np.zeros_like(mask)
    for i in range(1, n + 1):
        rr, cc = np.nonzero(lab == i)
        if len(cc) > 0 and (cc.max() - cc.min()) > 1200:
            keep[lab == i] = True
    return keep


def read_curve(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    mask = _curve_mask(image, which, r0, r1, c0, c1)
    return [
        C.sample_solid_curve_y(mask, xa, ya, float(x), r0, r1, note_prefix=f"FigS2C {which}: ")
        for x in at
    ]
