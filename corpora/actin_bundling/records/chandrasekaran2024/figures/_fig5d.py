"""Shared calibrated grouped-bar reader for Fig. 5D Percent Ring abundance vs T_2piR and [VASP-tet] (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

VASP_CONCS = [0.2, 0.4, 0.8, 1.0, 2.0]
T_2PIR_VALS = [10.0, 50.0, 150.0, 200.0, 450.0, 600.0]


def read_fig5d(image: np.ndarray, group_idx: int) -> list[dict]:
    """Read one [VASP-tet] group (0=0.2, 1=0.4, 2=0.8, 3=1.0, 4=2.0 uM) across all 6 T_2piR bars."""
    info = pc.find_plot_box(image)
    y_axis = pc.fit_axis(info["y_ticks"], [100.0, 80.0, 60.0, 40.0, 20.0, 0.0])
    x_ticks = info["x_ticks"]
    if len(x_ticks) != 30:
        raise ValueError(f"Expected 30 x-ticks in Fig. 5D, got {len(x_ticks)}")
    c0, r0, c1, r1 = info["box"]
    row_0 = y_axis.to_pixel(0.0)
    conc = VASP_CONCS[group_idx]

    readings = []
    for ti, tv in enumerate(T_2PIR_VALS):
        xt = x_ticks[group_idx * 6 + ti]
        cx = int(round(xt))
        strip = image[r0 + 2:r1 - 2, cx - 16:cx + 16].astype(np.float64)
        colored = (strip.max(axis=2) - strip.min(axis=2)) > 35.0
        rows = np.flatnonzero(colored.mean(axis=1) > 0.5) + (r0 + 2)
        if rows.size >= 5:
            # Find the center of the thin (~4 px) black top border immediately above the colored fill
            r_fill_top = int(rows.min())
            gray_above = image[max(r0, r_fill_top - 8):r_fill_top + 1, cx - 12:cx + 12].astype(np.float64).mean(axis=(1, 2))
            dark_b = np.flatnonzero(gray_above < 100.0)
            if dark_b.size:
                row_top = float(max(r0, r_fill_top - 8) + 0.5 * (dark_b[0] + dark_b[-1]) + 0.5)
            else:
                row_top = float(r_fill_top)
            note = f"Fig. 5D Ring abundance at [VASP-tet]={conc} uM, T_2piR={int(tv)} s (ry={y_axis.residual_px:.2f}px)"
            readings.append(pc.make_reading(tv, row_top, y_axis, xt, half_px=2.0, note=note, clamp_y=(0.0, 100.0)))
        else:
            # Zero-height bar sits cleanly at y = 0%
            note = f"Fig. 5D Ring abundance = 0% at [VASP-tet]={conc} uM, T_2piR={int(tv)} s (ry={y_axis.residual_px:.2f}px)"
            readings.append(pc.make_reading(tv, row_0, y_axis, xt, half_px=2.0, note=note, clamp_y=(0.0, 100.0)))
    return readings
