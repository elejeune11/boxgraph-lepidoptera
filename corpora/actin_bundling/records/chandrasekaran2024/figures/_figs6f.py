"""Shared calibrated curve reader for Supplementary Fig. 6F normalized radial density profiles (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

# Exact RGB colors for the distinguishable curves in Supp. Fig. 6F
SERIES_SPEC = {
    "nuc0": ((140.0, 200.0, 75.0), 32.0, "0 nM (green peripheral peak)"),
    "nuc75": ((196.0, 108.0, 42.0), 28.0, "75 nM (orange mid-droplet peak)"),
    "nuc150": ((156.0, 36.0, 45.0), 26.0, "150 nM (solid dark-maroon core peak)"),
}


def read_figs6f(image: np.ndarray, key: str, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], max_residual_px=4.5)
    y_axis = pc.fit_axis(info["y_ticks"], [30.0, 25.0, 20.0, 15.0, 10.0, 5.0, 0.0], max_residual_px=4.0)
    c0, r0, c1, r1 = info["box"]

    rgb, tol, desc = SERIES_SPEC[key]
    dist = pc.rgb_dist(image, rgb)
    mask = dist < tol
    mask[:r0 + 5, :] = False
    mask[r1 - 5:, :] = False
    mask[:, :c0 + 5] = False
    mask[:, c1 - 5:] = False
    # Exclude upper-right legend box
    mask[r0:r0 + 520, c0 + 360:c1] = False

    return pc.sample_curve_columns(
        mask, x_axis, y_axis, at, search_half_w=5, min_px_per_col=4,
        clamp_y=(0.0, 33.0), note_prefix=f"Supp. Fig. 6F [Nuc]={desc}",
    )
