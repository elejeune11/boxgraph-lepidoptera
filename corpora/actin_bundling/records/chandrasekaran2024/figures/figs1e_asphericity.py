"""Calibrated bar reader for Supplementary Fig. 1E LAMMPS droplet asphericity (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image, dark_thresh=160.0, min_tick_len=3, max_tick_w=30)
    y_axis = pc.fit_axis(info["y_ticks"], [1.2, 1.1, 1.0, 0.9, 0.8], max_residual_px=3.5)
    c0, r0, c1, r1 = info["box"]
    if len(info["x_ticks"]) != 2:
        raise ValueError(f"Expected 2 x-ticks in Supp. Fig. 1E, got {len(info['x_ticks'])}")

    readings = []
    for nf, xt in zip([6.0, 20.0], info["x_ticks"]):
        cx = int(round(xt))
        strip = image[r0 + 5:r1 - 5, cx - 100:cx + 100].astype(int)
        coral = (strip[:, :, 0] > 180) & (strip[:, :, 1] < 145) & (strip[:, :, 2] < 145)
        rows = np.flatnonzero(coral.mean(axis=1) > 0.5) + (r0 + 5)
        if not rows.size:
            raise ValueError(f"Could not locate coral bar at N={nf:g} in Supp. Fig. 1E")
        row_top = float(rows.min())
        note = f"Supp. Fig. 1E LAMMPS asphericity at N={int(nf)} filaments (ry={y_axis.residual_px:.2f}px)"
        readings.append(pc.make_reading(nf, row_top, y_axis, xt, half_px=2.0, note=note, clamp_y=(0.8, 1.2)))
    return readings
