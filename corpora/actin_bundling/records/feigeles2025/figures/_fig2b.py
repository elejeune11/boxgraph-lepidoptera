"""Shared calibrated reader for Fig. 2B (droplet aspect ratio vs time: blue dots + red exponential fit)."""

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
        C.find_ticks(dark, hb, (int(vl0[2]) + 15, int(vr0[1]) - 15), "x", side=-1, band=(8, 35)), 7
    )
    yt = C.match_ticks_regular(
        C.find_ticks(dark, vr0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=-1, band=(8, 35)), 4
    )
    xa = C.verify_axis(
        "linear", [(vl0[0], 0.0)] + list(zip(xt, [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0])) + [(vr0[0], 80.0)], 3.0
    )
    ya = C.verify_axis(
        "linear", [(ht[0], 1.1)] + list(zip(yt, [1.0, 0.9, 0.8, 0.7])) + [(hb[0], 0.6)], 3.0
    )
    return xa, ya, ht, hb, vl0, vr0


def read_data(image: np.ndarray) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    # Allow slight left-spine straddle for t=0 dot
    inb[int(ht[2]) + 6 : int(hb[1]) - 6, int(vl0[2]) - 15 : int(vr0[1]) - 6] = True
    # Exclude top-right inset cartoon (x > 30 s, y > 0.90)
    inb[: int(ya.to_pixel(0.90)), int(xa.to_pixel(30.0)) :] = False

    blue = inb & (B - R > 45) & (G < 110) & (B > 90)
    lab, n = ndimage.label(blue)
    areas = sorted([int((lab == i).sum()) for i in range(1, n + 1) if (lab == i).sum() > 200])
    # Single-dot radius from the median of isolated dots (excluding smaller stray speck / clipped dot)
    iso = [a for a in areas if 1900 <= a <= 3000]
    radius_px = math.sqrt(float(np.median(iso)) / math.pi)
    centers = C.split_disks_xor(
        blue, radius_px=radius_px, min_new_frac=0.14, min_sep_frac=0.20, min_xor_gain_frac=0.05, max_per_cluster=20
    )

    out = []
    for c, r, mult in centers:
        half_px = (2.5 if mult == 1 else 4.5) + ya.residual_px
        out.append(
            C.make_marker_reading(
                c,
                r,
                xa,
                ya,
                half_px,
                f"Fig2B blue circle (cluster={mult}), r={radius_px:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    out.sort(key=lambda d: d["x"])
    return out


def read_fit(image: np.ndarray, at: list[float]) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    r0, r1 = int(ht[2]) + 6, int(hb[1]) - 6
    c0, c1 = int(vl0[2]) + 6, int(vr0[1]) - 6
    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    red = inb & (R > 160) & (R - G > 70) & (R - B > 70)
    return [
        C.sample_curve_local_fit(red, xa, ya, xv, r0, r1, window_px=75, deg=2, note_prefix="Fig2B red exp fit: ")
        for xv in at
    ]
