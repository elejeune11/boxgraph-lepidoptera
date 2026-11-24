"""Calibrated reader for the three horizontal reference lines (C_total, C_total/2, C_total/3) in Supp. Fig. 13C (1200 dpi)."""

from __future__ import annotations
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [1.0, 2.0, 3.0, 4.0, 5.0, 5.56])
    y_axis = pc.fit_axis(info["y_ticks"], [80.0, 60.0, 40.0, 20.0, 0.0])
    c0, r0, c1, r1 = info["box"]

    gray = image[r0 + 10:r1 - 10, c0 + 10:c1 - 10].astype(np.float64).mean(axis=2)
    grey_frac = ((gray > 130.0) & (gray < 238.0)).mean(axis=1)
    hits = np.flatnonzero(grey_frac > 0.65) + (r0 + 10)
    if not hits.size:
        raise ValueError("Could not find grey reference lines in Supp. Fig. 13C")
    splits = np.flatnonzero(np.diff(hits) > 3)
    starts = np.concatenate(([0], splits + 1))
    ends = np.concatenate((splits, [hits.size - 1]))
    lines = [(0.5 * (float(hits[a]) + float(hits[b])) + 0.5, 0.5 * (float(hits[b]) - float(hits[a]) + 1.0))
             for a, b in zip(starts, ends)]
    if len(lines) != 3:
        raise ValueError(f"Expected 3 horizontal grey reference lines in Supp. Fig. 13C, got {len(lines)}")

    pcts = [100.0, 50.0, 33.34]
    labels = ["C_total (100%)", "C_total/2 (50%)", "C_total/3 (33.34%)"]
    col_mark = x_axis.to_pixel(3.5)
    readings = []
    for pct, lbl, (row_c, half_r) in zip(pcts, labels, lines):
        note = f"Supp. Fig. 13C horizontal reference line {lbl} (ry={y_axis.residual_px:.2f}px)"
        readings.append(pc.make_reading(pct, row_c, y_axis, col_mark, half_px=half_r, note=note))
    return readings
