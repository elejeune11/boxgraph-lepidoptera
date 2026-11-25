"""Calibrated reader for Main Fig. 4B (measured bisected daughter droplet lengths L_B vs L_A)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _calib import (
    detect_template_mp,
    find_blobs,
    find_box_spines,
    find_ticks,
    make_point_reading,
    match_ticks_regular,
    verify_axis,
)


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    dark, gray, hs, vs = find_box_spines(image, thresh=100, row_frac=0.45, col_frac=0.45)
    hs = [s for s in hs if (s[2] - s[1]) < 25]
    vs = [s for s in vs if (s[2] - s[1]) < 25]
    top, bot = hs[0], hs[-1]
    left, right = vs[0], vs[-1]

    # 4 interior x-ticks at 4, 5, 6, 7 um; left = 3, right = 8.
    xt = match_ticks_regular(find_ticks(dark, top, (left[2] + 10, right[1] - 10), "x", side=+1, band=(5, 25), frac=0.5), 4)
    x_marks = [(left[0], 3.0)] + [(xt[i], 4.0 + i) for i in range(4)] + [(right[0], 8.0)]
    x_axis = verify_axis("linear", x_marks, max_residual_px=2.0)

    # 4 interior y-ticks at 7, 6, 5, 4 um; top = 8, bot = 3.
    yt = match_ticks_regular(find_ticks(dark, left, (top[2] + 10, bot[1] - 10), "y", side=+1, band=(5, 25), frac=0.5), 4)
    y_marks = [(top[0], 8.0)] + [(yt[i], 7.0 - i) for i in range(4)] + [(bot[0], 3.0)]
    y_axis = verify_axis("linear", y_marks, max_residual_px=2.0)

    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    H, W = image.shape[:2]

    # Pure black dots (excludes blue y=x line where B > 100 and B - R > 30)
    roi = np.zeros((H, W), dtype=bool)
    roi[top[2] + 2 : bot[1] - 1, left[2] + 2 : right[1] - 1] = True
    # Exclude inset schematic of bisected droplet + L_A, L_B labels (top-left)
    roi[: int(y_axis.to_pixel(6.05)), : int(x_axis.to_pixel(5.45))] = False
    blk = (R < 50) & (G < 50) & (B < 50) & roi

    # Remove narrow axis ticks (9 px wide) via fast iterative 3x3 erosion/dilation, then restore dot edges
    blk_dots = ndimage.binary_dilation(ndimage.binary_erosion(blk, iterations=8), iterations=9) & blk

    # Measure exact isolated dot radius directly from single circles (area 1800..2100 px^2)
    iso = find_blobs(blk_dots, min_area=1750, max_area=2150)
    rad = float(np.sqrt(np.median([b["area"] for b in iso]) / np.pi))
    rad_i = int(round(rad))
    yy, xx = np.ogrid[-rad_i : rad_i + 1, -rad_i : rad_i + 1]
    tpl = ((xx * xx + yy * yy) <= rad * rad).astype(float)

    # Note: one dot sits right near bottom-left corner (3.04, 3.04) partially clipped by bottom spine,
    # and one near top-right (7.86, 7.95) partially clipped by top spine, so allow min_full_frac=0.55.
    found = detect_template_mp(
        blk_dots,
        tpl,
        min_full_frac=0.55,
        min_res_frac=0.38,
        min_dist=int(round(0.75 * rad)),
        refine_r=5,
        max_hits=100,
    )
    found.sort(key=lambda p: (p[1], -p[0]))
    cal_note = f"cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px"
    return [
        make_point_reading(c, r, x_axis, y_axis, half_px_y=4.0, note=f"dot R={rad:.1f}px, {cal_note}")
        for r, c in found
    ]
