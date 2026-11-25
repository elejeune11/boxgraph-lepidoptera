"""Shared calibrated reader for Fig. 3C (coalescing condensate length on bundle vs time)."""

from __future__ import annotations

import math
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
    xt = C.match_ticks_regular(
        C.find_ticks(dark, hb, (int(vl0[2]) + 15, int(vr0[1]) - 15), "x", side=-1, band=(8, 35)), 8
    )
    yt = C.match_ticks_regular(
        C.find_ticks(dark, vr0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=-1, band=(8, 35)), 6
    )
    # Right spine is at ~82 s (tick for 80 s sits inside right spine)
    xa = C.verify_axis(
        "linear", [(vl0[0], 0.0)] + list(zip(xt, [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0])), 3.0
    )
    ya = C.verify_axis(
        "linear", [(ht[0], 9.0)] + list(zip(yt, [8.5, 8.0, 7.5, 7.0, 6.5, 6.0])) + [(hb[0], 5.5)], 3.0
    )
    return xa, ya, ht, hb, vl0, vr0


def read_data(image: np.ndarray) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[int(ht[2]) + 32 : int(hb[1]) - 32, int(vl0[2]) + 6 : int(vr0[1]) - 4] = True
    # Exclude top-right cartoon (x > 30, y > 7.5) and text "tau = 14.4 s" (x in [28, 65], y in [6.3, 7.1])
    inb[: int(ya.to_pixel(7.5)), int(xa.to_pixel(30.0)) :] = False
    inb[int(ya.to_pixel(7.1)) : int(ya.to_pixel(6.3)), int(xa.to_pixel(28.0)) : int(xa.to_pixel(65.0))] = False

    blk = inb & (R < 60) & (G < 60) & (B < 60)
    # Exclude interior left/right tick marks away from the trajectory endpoints (y=8.5 tick sits below t=0 dot at 8.75)
    row_grid = np.arange(image.shape[0])[:, None]
    blk[:, : int(vl0[2]) + 35] &= row_grid < int(ya.to_pixel(8.62))
    blk[:, int(vr0[1]) - 35 :] &= row_grid > int(ya.to_pixel(6.0))

    lab, n = ndimage.label(blk)
    areas = sorted([int((lab == i).sum()) for i in range(1, n + 1) if (lab == i).sum() > 200])
    iso = [a for a in areas if 400 <= a <= 700]
    radius_px = math.sqrt(float(np.median(iso)) / math.pi)
    centers = C.split_disks_xor(
        blk, radius_px=radius_px, min_new_frac=0.14, min_sep_frac=0.20, min_xor_gain_frac=0.05, max_per_cluster=25
    )

    out = []
    for c, r, mult in centers:
        half_px = (2.0 if mult == 1 else 3.5) + ya.residual_px
        out.append(
            C.make_marker_reading(
                c,
                r,
                xa,
                ya,
                half_px,
                f"Fig3C black circle (cluster={mult}), r={radius_px:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    out.sort(key=lambda d: d["x"])
    return out


def read_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    r0, r1 = int(ht[2]) + 4, int(hb[1]) - 4
    c0, c1 = int(vl0[2]) + 4, int(vr0[1]) - 4
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    red = inb & (R > 160) & (R - G > 70) & (R - B > 70)
    return [
        C.sample_curve_local_fit(red, xa, ya, xv, r0, r1, window_px=75, deg=2, note_prefix="Fig3C red exp fit: ")
        for xv in at
    ]
