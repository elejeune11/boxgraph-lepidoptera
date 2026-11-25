"""Calibrated reader for Fig. S4 (contact angle vs droplet radius R_max, n = 28)."""

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
        C.find_ticks(dark, vl0, (int(ht[2]) + 15, int(hb[1]) - 15), "y", side=+1, band=(8, 35)), 4
    )
    xa = C.verify_axis("linear", [(vl0[0], 2.0)] + list(zip(xt, [3.0, 4.0, 5.0, 6.0])) + [(vr0[0], 7.0)], 3.0)
    ya = C.verify_axis("linear", [(ht[0], 60.0)] + list(zip(yt, [50.0, 40.0, 30.0, 20.0])) + [(hb[0], 10.0)], 3.0)

    inb = np.zeros(image.shape[:2], bool)
    inb[int(ht[2]) + 25 : int(hb[1]) - 25, int(vl0[2]) + 25 : int(vr0[1]) - 25] = True
    dots = inb & (gray < 100)
    lab, n = ndimage.label(dots)
    areas = [int((lab == i).sum()) for i in range(1, n + 1) if (lab == i).sum() > 500]
    radius_px = math.sqrt(float(np.median(areas)) / math.pi)
    centers = C.split_disks_xor(dots, radius_px=radius_px, min_new_frac=0.12, min_sep_frac=0.18, min_xor_gain_frac=0.05)

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
                f"FigS4 black circle (cluster={mult}), r={radius_px:.1f}px, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    out.sort(key=lambda d: d["x"])
    return out
