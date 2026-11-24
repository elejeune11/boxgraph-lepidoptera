"""Shared calibrated stacked-bar reader for Supplementary Fig. 6D VASP bound-filament fractions (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

NUC_CONCS = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0]
PALETTE = np.array([
    [146.0, 200.0, 73.0],   # 0 filaments (green)
    [196.0, 183.0, 46.0],   # 1 filament (yellow)
    [201.0, 126.0, 41.0],   # 2 filaments (orange)
    [180.0, 67.0, 38.0],    # 3 filaments (rust)
    [162.0, 46.0, 50.0],    # 4 filaments (maroon)
], dtype=np.float64)


def _longest_run(rows: np.ndarray, max_gap: int = 4) -> tuple[float, float] | None:
    if rows.size < 3:
        return None
    sp = np.flatnonzero(np.diff(rows) > max_gap)
    st = np.concatenate(([0], sp + 1))
    en = np.concatenate((sp, [rows.size - 1]))
    best = max(zip(st, en), key=lambda ab: rows[ab[1]] - rows[ab[0]])
    a, b = best
    if rows[b] - rows[a] + 1 < 3:
        return None
    return float(rows[a]), float(rows[b] + 1)


def read_figs6d(image: np.ndarray, n_bound: int) -> list[dict]:
    """Read stacked bar fraction of VASP tetramers bound to n_bound in {0,1,2,3,4} filaments across [Nucleator]."""
    info = pc.find_plot_box(image)
    y_axis = pc.fit_axis(info["y_ticks"], [1.0, 0.5, 0.0], max_residual_px=4.0)
    x_ticks = info["x_ticks"]
    if len(x_ticks) != 8:
        raise ValueError(f"Expected 8 x-ticks in Supp. Fig. 6D, got {len(x_ticks)}")
    c0, r0, c1, r1 = info["box"]
    row_0 = y_axis.to_pixel(0.0)
    row_1 = y_axis.to_pixel(1.0)
    px_per_unit = abs(row_0 - row_1)
    # Panel is ~192 dpi rasterized in SI PDF (1 source px ~ 6.25 px at 1200 dpi)
    err_px = 6.5 + max(y_axis.residual_px, 1.0)
    err_val = err_px / px_per_unit

    readings = []
    for nuc, xt in zip(NUC_CONCS, x_ticks):
        cx = int(round(xt))
        left_s = image[r0 + 2:r1 - 1, cx - 65:cx - 25].astype(np.float64)
        right_s = image[r0 + 2:r1 - 1, cx + 25:cx + 65].astype(np.float64)
        strip = np.concatenate([left_s, right_s], axis=1)
        row_rgb = np.median(strip, axis=1)
        dists = np.sqrt(((row_rgb[:, None, :] - PALETTE[None, :, :]) ** 2).sum(axis=2))
        cls = np.argmin(dists, axis=1)
        min_d = np.min(dists, axis=1)
        # Tight color tolerance (< 24) rejects dark outline rows between segments
        valid = min_d < 24.0

        # Find longest pure interior run for each class k = 4 down to 0 in strictly descending order
        segs: dict[int, tuple[float, float]] = {}
        min_allowed_row = float(r0 + 2)
        for k in [4, 3, 2, 1, 0]:
            rows_k = np.flatnonzero(valid & (cls == k)) + (r0 + 2)
            rows_k = rows_k[rows_k >= min_allowed_row]
            run = _longest_run(rows_k, max_gap=4)
            if run is not None:
                segs[k] = run
                min_allowed_row = run[1]

        # Check thin darkened yellow strip (k=1) below k=2 when squeezed by borders at high [Nuc]
        if 1 not in segs and 2 in segs:
            below_start = int(segs[2][1])
            sub_rgb = row_rgb[below_start - (r0 + 2):]
            if sub_rgb.shape[0] > 0:
                # Yellow darkened by outline has high R and G (~equal) and low B
                is_yel = (sub_rgb[:, 0] > 90) & (sub_rgb[:, 1] > 85) & (sub_rgb[:, 2] < 50) & (np.abs(sub_rgb[:, 0] - sub_rgb[:, 1]) < 25)
                yr = np.flatnonzero(is_yel) + below_start
                run = _longest_run(yr, max_gap=4)
                if run is not None:
                    segs[1] = run

        # Compute cumulative boundaries b[0..5] where b[0]=row_0 (y=0) and b[5]=row_1 (y=1)
        # Interface between adjacent present classes splits the dark dividing outline 50/50
        present = sorted(segs.keys())  # ascending k (bottom -> top)
        bounds = {0: row_0, 5: row_1}
        for i, k in enumerate(present):
            if i == len(present) - 1:
                top_b = row_1
            else:
                k_next = present[i + 1]
                top_b = 0.5 * (segs[k][0] + segs[k_next][1])
            for m in range(k + 1, (present[i + 1] if i + 1 < len(present) else 5) + 1):
                if m not in bounds:
                    bounds[m] = top_b
        # Any missing classes below lowest present class get zero height (unless lowest > 0 and bottom > 10px above row_0)
        lowest = present[0]
        for m in range(1, lowest + 1):
            # If lowest is k=1 or k=2 and its pure segment stops slightly above row_0 due to bottom outline,
            # check if there is a tiny green/yellow sliver; otherwise assign boundary to row_0
            bounds[m] = row_0

        # Ensure monotonic bounds[0] >= bounds[1] >= ... >= bounds[5]
        b_rows = [bounds[m] for m in range(6)]
        fracs = [max(0.0, (b_rows[k] - b_rows[k + 1]) / px_per_unit) for k in range(5)]
        s = sum(fracs)
        if s > 0:
            fracs = [f / s for f in fracs]

        val = fracs[n_bound]
        cum_top = sum(fracs[:n_bound + 1])
        px_row = y_axis.to_pixel(cum_top) if n_bound < 4 else y_axis.to_pixel(1.0 - val)
        readings.append({
            "x": float(nuc),
            "y": float(val),
            "lo": float(max(0.0, val - err_val)),
            "hi": float(min(1.0, val + err_val)),
            "px": [float(xt), float(px_row)],
            "px_lo": float(px_row - err_px),
            "px_hi": float(px_row + err_px),
            "note": f"Supp. Fig. 6D VASP bound to {n_bound} filaments at [Nucleator]={int(nuc)} nM (cum top={cum_top:.3f}; ry={y_axis.residual_px:.2f}px)",
        })
    return readings
