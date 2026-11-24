"""Shared calibrated reader for Main Figs. 3e and 3f dashed power-law trendlines."""

from __future__ import annotations

import numpy as np
from scipy import ndimage
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_l_spines(image, row_frac=0.4, col_frac=0.4)
    hb = hl[-1]
    vl0 = vl[0]
    xt = C.find_outward_ticks(dark, hb, (int(vl0[1]) - 10, image.shape[1]), "x", side=+1, band=(6, 25))
    yt = C.find_outward_ticks(dark, vl0, (0, int(hb[2]) + 10), "y", side=-1, band=(6, 25))
    xt = C.match_ticks_regular(xt, 7, tol=6.0)
    yt = C.match_ticks_regular(yt, 10, tol=6.0)
    xa = C.verify_axis("linear", list(zip(xt, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])), max_residual_px=5.0)
    ya = C.verify_axis("linear", list(zip(yt, [1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0])), max_residual_px=5.0)
    return xa, ya, int(yt[0]) - 20, int(hb[1]) - 3


def read_power_fit(image: np.ndarray, at: list[float], color: str, label: str) -> list[dict]:
    xa, ya, row_min, row_max = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)

    if color == "blue":
        # Dark blue dashed fit line RGB ~(0, 112, 192), distinct from pale sky-blue scatter dots (~155, 206, 245)
        target = np.array([0.0, 112.0, 192.0])
        dist = np.sqrt((r - target[0]) ** 2 + (g - target[1]) ** 2 + (b - target[2]) ** 2)
        mask = (dist < 65.0) & (r < 70.0) & (b > 130.0)
    elif color == "red":
        # Pure red dashed fit line (~(192,0,0) in Fig. 3e, ~(255,0,0) in Fig. 3f).
        # Strip circular red scatter dots (~32px diameter in Fig. 3f) via morphological opening, leaving only thin dashes (~14px wide).
        raw_red = (r > 155.0) & (g < 35.0) & (b < 35.0)
        yy, xx = np.ogrid[-10:11, -10:11]
        disk = (xx * xx + yy * yy) <= 10 * 10
        dots = ndimage.binary_dilation(ndimage.binary_opening(raw_red, structure=disk), iterations=1)
        mask = raw_red & ~dots
    else:
        raise ValueError(color)

    # Mask upper-left legend area and right micrograph insets (x > 5.45)
    mask[: int(ya.to_pixel(1.0)), : int(xa.to_pixel(2.8))] = False
    mask[:, int(xa.to_pixel(5.45)) :] = False

    return [
        C.sample_dashed_curve_y(mask, xa, ya, x, row_min, row_max, max_search_dx=125, note_prefix=f"{label}: ")
        for x in at
    ]
