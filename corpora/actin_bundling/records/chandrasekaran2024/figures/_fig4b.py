"""Shared calibrated curve reader for Fig. 4B.i and 4B.ii pairwise filament-end distance PDFs (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

PALETTE_I = [
    (135.0, 186.0, 78.0),   # pi/4
    (163.0, 187.0, 50.0),   # pi/3
    (199.0, 187.0, 44.0),   # pi/2
    (214.0, 168.0, 38.0),   # 2pi/3
]
PALETTE_II = [
    (203.0, 126.0, 39.0),   # pi
    (194.0, 86.0, 40.0),    # 4pi/3
    (183.0, 65.0, 37.0),    # 5pi/3
    (173.0, 52.0, 40.0),    # 2pi
]


def read_fig4b(image: np.ndarray, panel: str, curve_idx: int, label: str, at: list[float]) -> list[dict]:
    """Read one curve from Fig. 4B.i or 4B.ii via nearest-prototype RGB classification."""
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [0.0, 0.5, 1.0, 1.5, 2.0])
    y_axis = pc.fit_axis(info["y_ticks"], [2.0, 1.5, 1.0, 0.5, 0.0])
    c0, r0, c1, r1 = info["box"]

    palette = PALETTE_I if panel == "i" else PALETTE_II
    arr = image.astype(np.float64)
    dists = np.stack([np.sqrt(np.sum((arr - np.array(rgb)) ** 2, axis=2)) for rgb in palette], axis=0)
    best_idx = np.argmin(dists, axis=0)
    min_d = np.min(dists, axis=0)

    mask = (best_idx == curve_idx) & (min_d < 22.0)
    mask[:r0 + 8, :] = False
    mask[r1 - 8:, :] = False
    mask[:, :c0 + 8] = False
    mask[:, c1 - 8:] = False
    # Exclude legend swatches in top-left corner of box
    mask[r0:r0 + 740, c0:c0 + 620] = False

    return pc.sample_curve_columns(
        mask, x_axis, y_axis, at, search_half_w=4, min_px_per_col=5,
        clamp_y=(0.0, 2.25), note_prefix=f"Fig. 4B.{panel} L_fil^max={label}",
    )
