"""Shared calibrated reader for Fig. 7F 2D droplet Bending (green) and Surface (brown) energy curves (1200 dpi)."""

from __future__ import annotations
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc


def read_fig7f(image: np.ndarray, curve: str, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [1.0, 3.0, 5.0, 7.0, 9.0, 11.0])
    y_axis = pc.fit_axis(info["y_ticks"], [20.0, 10.0, 0.0])
    c0, r0, c1, r1 = info["box"]

    arr = image.astype(int)
    if curve == "bending":
        mask = (arr[:, :, 1] - arr[:, :, 0] > 40) & (arr[:, :, 1] - arr[:, :, 2] > 20) & (arr[:, :, 1] > 80)
        label = "Fig. 7F Bending energy (green)"
    elif curve == "surface":
        mask = (
            (arr[:, :, 0] - arr[:, :, 1] > 12)
            & (arr[:, :, 0] - arr[:, :, 2] > 12)
            & (arr[:, :, 0] > 70)
            & (arr[:, :, 0] < 140)
            & (arr[:, :, 1] < 95)
        )
        label = "Fig. 7F Surface energy (brown)"
    else:
        raise ValueError(f"Unknown curve {curve!r}")

    # Restrict strictly to plot box and mask out legend swatches in upper middle
    mask[:r0 + 5, :] = False
    mask[r1 - 5:, :] = False
    mask[:, :c0 + 5] = False
    mask[:, c1 - 5:] = False
    mask[r0:r0 + 580, c0 + 380:c0 + 1350] = False

    return pc.sample_curve_columns(mask, x_axis, y_axis, at, search_half_w=2, note_prefix=label)
