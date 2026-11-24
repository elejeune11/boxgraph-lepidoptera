"""Shared calibrated readers for PCA cumulative explained variance (S4A, S9A, S11Bi, S12A)
and k-means Silhouette bar charts (S4C, S9C, S11Biii, S12C)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc


def read_pca_cumvar(image: np.ndarray, n_pcs: int, panel_label: str) -> list[dict]:
    """Read cumulative explained variance ratio after PC1..PC_{n_pcs} from horizontal stacked bar."""
    info = pc.find_plot_box(image)
    c0, r0, c1, r1 = info["box"]
    x_ticks = [t for t in info["x_ticks"] if t < info["spines"]["right"] - 15]
    x_vals = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    x_axis = pc.fit_axis(x_ticks, x_vals, max_residual_px=3.5)

    # Locate vertical extent of the dark-maroon PC1 bar at x = 0.1
    cx01 = int(round(x_axis.to_pixel(0.1)))
    col01 = image[r0 + 5:r1 - 5, cx01].astype(np.float64)
    bar_rows = np.flatnonzero((col01[:, 0] > 60.0) & (col01[:, 1] < 45.0)) + (r0 + 5)
    if not bar_rows.size:
        raise ValueError(f"Could not find PC1 bar in {panel_label}")
    r_lo, r_hi = int(bar_rows.min()), int(bar_rows.max())
    row_mid = 0.5 * (r_lo + r_hi) + 0.5

    # Clean bottom strip inside bar (free of centered PC labels and upper arrows)
    rs0 = int(r_hi - 0.18 * (r_hi - r_lo))
    rs1 = int(r_hi - 0.04 * (r_hi - r_lo))
    strip = image[rs0:rs1, c0:c1].astype(np.float64).mean(axis=0)
    diffs = np.linalg.norm(np.diff(strip, axis=0), axis=1)

    # Find sharp color transitions strictly within (0.05, 1.02]
    hits = np.flatnonzero(diffs > 6.0)
    if not hits.size:
        raise ValueError(f"No color transitions found in {panel_label}")
    splits = np.flatnonzero(np.diff(hits) > 5)
    gs = np.concatenate(([0], splits + 1))
    ge = np.concatenate((splits, [hits.size - 1]))
    trans = []
    for a, b in zip(gs, ge):
        idx = hits[a:b + 1]
        pk = idx[int(np.argmax(diffs[idx]))]
        col_abs = float(c0 + pk + 1.0)
        val = x_axis.to_value(col_abs)
        if 0.05 < val <= 1.015:
            trans.append((val, col_abs))
    if len(trans) < n_pcs:
        raise ValueError(f"{panel_label}: expected {n_pcs} transitions, found {len(trans)} ({[round(t[0],4) for t in trans]})")
    # Keep the first n_pcs transitions (PC1..PC_{n_pcs})
    trans = trans[:n_pcs]

    px_per_unit = abs(x_axis.to_pixel(1.0) - x_axis.to_pixel(0.0))
    err_val = (2.0 + max(x_axis.residual_px, 0.5)) / px_per_unit

    readings = []
    for k, (val, col_abs) in enumerate(trans, start=1):
        vc = min(1.0, max(0.0, val))
        lo = max(0.0, vc - err_val)
        hi = min(1.0, vc + err_val)
        readings.append({
            "x": float(k),
            "y": float(vc),
            "lo": float(lo),
            "hi": float(hi),
            "px": [float(col_abs), float(row_mid)],
            "px_lo": float(r_lo),
            "px_hi": float(r_hi),
            "note": f"{panel_label} cumulative explained variance through PC{k} (rx={x_axis.residual_px:.2f}px)",
        })
    return readings


def read_silhouette(image: np.ndarray, y_tick_vals: list[float], x_tick_vals: list[float],
                    panel_label: str) -> list[dict]:
    """Read all 9 k-means Silhouette bars (k = 2..10)."""
    info = pc.find_plot_box(image)
    c0, r0, c1, r1 = info["box"]
    x_ticks = [t for t in info["x_ticks"] if c0 + 20 < t < c1 - 20]
    x_axis = pc.fit_axis(x_ticks, x_tick_vals, max_residual_px=4.5)
    y_axis = pc.fit_axis(info["y_ticks"], y_tick_vals, max_residual_px=3.0)

    readings = []
    for k in range(2, 11):
        pc_x = x_axis.to_pixel(float(k))
        cx = int(round(pc_x))
        hw = max(12, int(round(0.15 * abs(x_axis.to_pixel(3.0) - x_axis.to_pixel(2.0)))))
        strip = image[r0 + 4:r1 - 2, max(c0, cx - hw):min(c1, cx + hw)].astype(np.float64)
        colored = (strip.max(axis=2) - strip.min(axis=2)) > 25.0
        rows = np.flatnonzero(colored.mean(axis=1) > 0.5) + (r0 + 4)
        if not rows.size:
            raise ValueError(f"{panel_label}: could not find bar for k={k}")
        row_top = float(rows.min())
        note = f"{panel_label} silhouette score at k={k} (rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(float(k), row_top, y_axis, pc_x, half_px=1.5, x_axis=x_axis, note=note, clamp_y=(0.0, 1.0))
        )
    return readings
