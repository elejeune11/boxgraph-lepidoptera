"""Fig. S2 (PDF p. 16 / SI-10): growth exponent alpha vs. normalized actin intensity I/I_max."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import build_axis, find_box, match_ticks_both_sides, reading_from_px, scan_ticks


def read(image, at):
    box = find_box(image)
    c0, r0, c1, r1 = box

    # 5 interior x-ticks at 0.2, 0.4, 0.6, 0.8, 1.0 plus box edges 0.0 and 1.2
    xt = match_ticks_both_sides(
        scan_ticks(image, box, "bottom_in", offset_min=8, offset_max=35, min_len=12),
        scan_ticks(image, box, "top_in", offset_min=8, offset_max=35, min_len=12),
    )
    x_pos = [c0] + [p for p, _ in xt] + [c1]
    x_vals = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    xa = build_axis("linear", x_pos, x_vals, max_residual_px=2.0)

    # 4 interior y-ticks at 0.70, 0.65, 0.60, 0.55 plus box edges 0.75 and 0.50
    yt = match_ticks_both_sides(
        scan_ticks(image, box, "left_in", offset_min=8, offset_max=35, min_len=12),
        scan_ticks(image, box, "right_in", offset_min=8, offset_max=35, min_len=12),
    )
    y_pos = [r0] + [p for p, _ in yt] + [r1]
    y_vals = [0.75, 0.70, 0.65, 0.60, 0.55, 0.50]
    ya = build_axis("linear", y_pos, y_vals, max_residual_px=2.0)

    # Isolate dark neutral ink inside the plot box (avoiding frame/ticks)
    rgb = image.astype(np.float64)
    gray = rgb.mean(axis=2)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    dark = (gray < 120.0) & (chroma < 35.0)
    mask = np.zeros_like(dark)
    mask[int(r0) + 40 : int(r1) - 40, int(c0) + 40 : int(c1) - 40] = dark[
        int(r0) + 40 : int(r1) - 40, int(c0) + 40 : int(c1) - 40
    ]

    # Identify open circular markers purely by ring geometry (square bbox ~60-110 px + enclosed white hole)
    lbl, n = ndimage.label(mask)
    slices = ndimage.find_objects(lbl)
    markers = []
    for i, sl in enumerate(slices, start=1):
        if sl is None:
            continue
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if not (55 <= h <= 130 and 55 <= w <= 130):
            continue
        if max(h, w) / min(h, w) > 1.15:
            continue
        sub = lbl[sl] == i
        area = int(sub.sum())
        filled = ndimage.binary_fill_holes(sub)
        hole_area = int(filled.sum()) - area
        if hole_area < 0.4 * area:
            continue
        rr, cc = np.nonzero(filled)
        row = float(rr.mean()) + sl[0].start + 0.5
        col = float(cc.mean()) + sl[1].start + 0.5
        markers.append((col, row, w, h, area))

    markers.sort(key=lambda t: t[0])
    if len(markers) != 5:
        raise RuntimeError(f"Fig. S2: expected 5 open circular markers, found {len(markers)}")

    out = []
    for col, row, w, h, area in markers:
        out.append(
            reading_from_px(
                xa,
                ya,
                col,
                row,
                6.0,
                note=f"open circle centroid ({w}x{h} px, ring {area} px, cal res ({xa.residual_px:.2f}, {ya.residual_px:.2f}) px)",
            )
        )
    return out
