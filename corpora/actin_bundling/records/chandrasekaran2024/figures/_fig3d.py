"""Shared calibrated stacked-bar reader for Fig. 3D VASP bound-filament fractions across all 6 k_bind panels (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

K_BIND_LABELS = ["10.0/s", "1.0/s", "0.1/s", "0.01/s", "0.001/s", "0.0001/s"]
K_UNBIND_VALS = [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0]
PALETTE = np.array([
    [135.0, 186.0, 78.0],   # 0 filaments (green)
    [191.0, 179.0, 43.0],   # 1 filament (yellow)
    [203.0, 126.0, 39.0],   # 2 filaments (orange)
    [183.0, 65.0, 37.0],    # 3 filaments (rust)
    [165.0, 44.0, 48.0],    # 4 filaments (maroon)
], dtype=np.float64)


def read_fig3d(image: np.ndarray, panel_idx: int, n_bound: int) -> list[dict]:
    """Read stacked fraction of VASP tetramers bound to n_bound in {0..4} filaments for one k_bind row (0=10..5=1e-4)."""
    boxes = pc.find_stacked_boxes(image, 6)
    b = boxes[panel_idx]
    c0, r0, c1, r1 = b["box"]
    y_axis = pc.fit_axis(b["y_ticks"], [1.0, 0.5, 0.0], max_residual_px=3.0)
    # Use bottom panel x_ticks (shared vertically across all 6 aligned subpanels)
    x_ticks = boxes[5]["x_ticks"]
    if len(x_ticks) != 6:
        raise ValueError(f"Expected 6 x-ticks in Fig. 3D, got {len(x_ticks)}")
    row_0 = y_axis.to_pixel(0.0)
    row_1 = y_axis.to_pixel(1.0)
    px_per_unit = abs(row_0 - row_1)
    err_val = (3.0 + max(y_axis.residual_px, 0.5)) / px_per_unit
    kb_lbl = K_BIND_LABELS[panel_idx]

    readings = []
    for ku, xt in zip(K_UNBIND_VALS, x_ticks):
        cx = int(round(xt))
        left_s = image[r0 + 2:r1 - 1, cx - 75:cx - 25].astype(np.float64)
        right_s = image[r0 + 2:r1 - 1, cx + 25:cx + 75].astype(np.float64)
        strip = np.concatenate([left_s, right_s], axis=1)
        row_rgb = np.median(strip, axis=1)
        dists = np.sqrt(((row_rgb[:, None, :] - PALETTE[None, :, :]) ** 2).sum(axis=2))
        cls = np.argmin(dists, axis=1)
        min_d = np.min(dists, axis=1)
        valid = min_d < 42.0

        fracs = []
        for k in range(5):
            rows_k = np.flatnonzero(valid & (cls == k)) + (r0 + 2)
            if rows_k.size >= 3:
                top_r = float(rows_k.min()) if k < 4 else row_1
                bot_r = float(rows_k.max() + 1) if k > 0 else row_0
                fracs.append(max(0.0, (bot_r - top_r) / px_per_unit))
            else:
                fracs.append(0.0)
        s = sum(fracs)
        if s > 0:
            fracs = [f / s for f in fracs]

        val = fracs[n_bound]
        cum_top = sum(fracs[:n_bound + 1])
        px_row = y_axis.to_pixel(cum_top) if n_bound < 4 else y_axis.to_pixel(1.0 - val)
        h_px = err_val * px_per_unit
        readings.append({
            "x": float(ku),
            "y": float(val),
            "lo": float(max(0.0, val - err_val)),
            "hi": float(min(1.0, val + err_val)),
            "px": [float(xt), float(px_row)],
            "px_lo": float(px_row - h_px),
            "px_hi": float(px_row + h_px),
            "note": f"Fig. 3D k_bind={kb_lbl}, k_unbind={ku:g}/s VASP bound to {n_bound} filaments (cum top={cum_top:.3f}; ry={y_axis.residual_px:.2f}px)",
        })
    return readings
