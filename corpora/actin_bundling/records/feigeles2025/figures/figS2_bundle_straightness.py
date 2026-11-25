"""Calibrated reader for Fig. S2 (bundle straightness ratio vs contour length L_c, n = 50)."""

from __future__ import annotations

import math
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hl, vl = C.find_box_spines(image, thresh=140.0)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.match_ticks_regular(
        C.find_ticks(dark, hb, (int(vl0[2]) + 15, int(vr0[1]) - 15), "x", side=-1, band=(8, 35)), 4
    )
    yt = C.match_ticks_regular(
        C.find_ticks(dark, vr0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=-1, band=(8, 35)), 4
    )
    xa = C.verify_axis("linear", [(vl0[0], 0.0)] + list(zip(xt, [10.0, 20.0, 30.0, 40.0])) + [(vr0[0], 50.0)], 3.0)
    ya = C.verify_axis("linear", [(ht[0], 3.0)] + list(zip(yt, [2.5, 2.0, 1.5, 1.0])) + [(hb[0], 0.5)], 3.0)

    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[int(ht[2]) + 25 : int(hb[1]) - 25, int(vl0[2]) + 25 : int(vr0[1]) - 25] = True
    # Exclude upper-left schematic inset (x < 35 um, y > 1.8)
    inb[: int(ya.to_pixel(1.8)), : int(xa.to_pixel(35.0))] = False

    blk = inb & (R < 60) & (G < 60) & (B < 60)
    lab, n = ndimage.label(blk)
    areas = sorted([int((lab == i).sum()) for i in range(1, n + 1) if (lab == i).sum() > 300])
    iso = [a for a in areas if 1500 <= a <= 2500]
    radius_px = math.sqrt(float(np.median(iso)) / math.pi)
    centers = C.split_disks_xor(
        blk, radius_px=radius_px, min_new_frac=0.12, min_sep_frac=0.18, min_xor_gain_frac=0.04, max_per_cluster=14
    )

    out = []
    for c, r, mult in centers:
        half_px = (2.0 if mult == 1 else 4.0) + ya.residual_px
        out.append(
            C.make_marker_reading(
                c,
                r,
                xa,
                ya,
                half_px,
                f"FigS2 black circle (cluster={mult}), r={radius_px:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    out.sort(key=lambda d: d["x"])
    return out
