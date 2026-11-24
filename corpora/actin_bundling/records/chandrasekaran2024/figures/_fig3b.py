"""Shared calibrated curve reader for distinguishable surface-area-fraction kinetics in Fig. 3B (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

K_BIND_LABELS = ["10.0/s", "1.0/s", "0.1/s", "0.01/s", "0.001/s", "0.0001/s"]
# Distinct solid mean-line colors (darker & more saturated than their semi-transparent SD bands)
COLOR_SPECS = {
    "0.0001": ((90.0, 158.0, 72.0), 38.0),    # bright green
    "1.0":    ((171.0, 99.0, 40.0), 16.0),    # burnt orange/brown (ring trajectory; distinct from 0.1/s=(163,122,43))
    "10.0":   ((174.0, 72.0, 37.0), 20.0),    # brick red (topmost drawn curve)
}


def read_fig3b(image: np.ndarray, panel_idx: int, ku_key: str, at: list[float]) -> list[dict]:
    """Read one cleanly separated surface-area-fraction trajectory from Fig. 3B."""
    boxes = pc.find_stacked_boxes(image, 6)
    b = boxes[panel_idx]
    c0, r0, c1, r1 = b["box"]
    x_axis = pc.fit_axis(boxes[5]["x_ticks"], [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0], max_residual_px=3.0)
    y_axis = pc.fit_axis(b["y_ticks"], [0.6, 0.3, 0.0], max_residual_px=3.0)

    rgb, tol = COLOR_SPECS[ku_key]
    dist = pc.rgb_dist(image, rgb)
    # Also require darker intensity than faint shaded SD bands
    lum = image.astype(np.float64).mean(axis=2)
    mask = (dist < tol) & (lum < 155.0)
    mask[:r0 + 3, :] = False
    mask[r1 - 2:, :] = False
    mask[:, :c0 + 3] = False
    mask[:, c1 - 3:] = False

    prefer = "top" if ku_key == "0.0001" else "median"
    return pc.sample_curve_columns(
        mask, x_axis, y_axis, at, search_half_w=6, prefer=prefer, min_px_per_col=4,
        clamp_y=(0.0, 0.6), note_prefix=f"Fig. 3B k_bind={K_BIND_LABELS[panel_idx]}, k_unbind={ku_key}/s",
    )
