"""Shared calibrated marker reader for Supplementary Fig. 1B radial density profiles (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy import ndimage
import _plot_common as pc
from boxgraph import digitize as D


def read_figs1b(image: np.ndarray, series: str) -> list[dict]:
    """Read circular markers of 'vasp' (red Droplet) or 'crowder' (black Crowder) in Supp. Fig. 1B."""
    info = pc.find_plot_box(image, dark_thresh=215.0, min_tick_len=3, max_tick_w=30)
    c0, r0, c1, r1 = info["box"]
    x_axis = pc.fit_axis(info["x_ticks"], [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0], max_residual_px=3.5)

    # Match detected y-ticks to regular [0.8, 0.7, ..., 0.0] grid robustly
    yt = info["y_ticks"]
    if len(yt) == 9:
        y_axis = pc.fit_axis(yt, [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0], max_residual_px=4.0)
    else:
        step = (yt[-1] - yt[0]) / 8.0
        pairs = [(t, round(0.8 - 0.1 * round((t - yt[0]) / step), 1)) for t in yt]
        y_axis = D.Axis("linear", pairs)

    sub = image.astype(int)
    if series == "vasp":
        mask = (sub[:, :, 0] > 180) & (sub[:, :, 1] < 95) & (sub[:, :, 2] < 95)
        label = "Supp. Fig. 1B Droplet (VASP) radial density"
    elif series == "crowder":
        mask = (sub.max(axis=2) < 85) & ((sub.max(axis=2) - sub.min(axis=2)) < 30)
        label = "Supp. Fig. 1B Crowder radial density"
    else:
        raise ValueError(f"Unknown series {series!r}")

    # Restrict strictly to plot box and exclude legend box + panel label 'B.'
    mask[:r0 + 8, :] = False
    mask[r1 - 8:, :] = False
    mask[:, :c0 + 8] = False
    mask[:, c1 - 8:] = False
    mask[r0 + 540:r0 + 1040, c0 + 720:c1 - 8] = False  # legend box
    mask[r0:r0 + 240, c0:c0 + 220] = False              # 'B.' label

    # Isolate filled circular markers (diameter ~28 px; connecting lines are ~6 px thick)
    dots = ndimage.binary_opening(mask, structure=np.ones((15, 15)))
    lab, n = ndimage.label(dots)
    pts = []
    for k in range(1, n + 1):
        rr, cc = np.nonzero(lab == k)
        if rr.size >= 150:
            r_c = float(rr.mean()) + 0.5
            c_c = float(cc.mean()) + 0.5
            pts.append((x_axis.to_value(c_c), c_c, r_c, len(rr)))
    pts.sort(key=lambda t: t[0])

    readings = []
    for xv, c_c, r_c, area in pts:
        note = f"{label} (marker {area}px; rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(xv, r_c, y_axis, c_c, half_px=2.0, x_axis=x_axis, note=note, clamp_y=(0.0, 0.85))
        )
    return readings
