"""Shared calibrated curve reader for Supplementary Fig. 6B filament count time series (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

NUC_CONCS = [0, 10, 25, 50, 75, 100, 125, 150]


def read_figs6b(image: np.ndarray, curve_idx: int, at: list[float]) -> list[dict]:
    """Read one curve (curve_idx 0..7 corresponding to [Nucleator]=0..150 nM, strictly ordered bottom->top)."""
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0], max_residual_px=4.0)
    y_axis = pc.fit_axis(info["y_ticks"], [1500.0, 1000.0, 500.0, 0.0], max_residual_px=4.0)
    c0, r0, c1, r1 = info["box"]

    arr = image.astype(int)
    # Solid mean lines have higher saturation and darker minimum channel than the faint shaded SD bands
    solid = ((arr.max(axis=2) - arr.min(axis=2)) > 75) & (arr.min(axis=2) < 110)
    solid[:r0 + 5, :] = False
    solid[r1 - 5:, :] = False
    solid[:, :c0 + 5] = False
    solid[:, c1 - 5:] = False
    solid[r0:r0 + 550, c0:c0 + 1150] = False  # legend box

    conc = NUC_CONCS[curve_idx]
    readings = []
    for tv in at:
        pc_x = x_axis.to_pixel(tv)
        col = int(round(pc_x))
        rs = np.flatnonzero(solid[:, col])
        if not rs.size:
            readings.append({"x": float(tv), "y": None, "note": f"no curves at t={tv:g} s"})
            continue
        sp = np.flatnonzero(np.diff(rs) > 10)
        st = np.concatenate(([0], sp + 1))
        en = np.concatenate((sp, [rs.size - 1]))
        runs = [(0.5 * (float(rs[a]) + float(rs[b])) + 0.5, 0.5 * (float(rs[b]) - float(rs[a]) + 1.0))
                for a, b in zip(st, en) if (rs[b] - rs[a] + 1) >= 6]
        # Reverse so index 0 is bottommost (0 nM) and index 7 is topmost (150 nM)
        runs.reverse()
        if len(runs) != 8:
            readings.append({"x": float(tv), "y": None, "note": f"expected 8 curves at t={tv:g} s, found {len(runs)}"})
            continue
        row_c, half_r = runs[curve_idx]
        # Local slope from +-6 px columns
        slope = 0.0
        if c0 + 10 < col - 6 and col + 6 < c1 - 10:
            r_l = np.flatnonzero(solid[:, col - 6])
            r_r = np.flatnonzero(solid[:, col + 6])
            rl_near = r_l[np.abs(r_l + 0.5 - row_c) < 35.0]
            rr_near = r_r[np.abs(r_r + 0.5 - row_c) < 35.0]
            if rl_near.size and rr_near.size:
                slope = (float(np.median(rr_near)) - float(np.median(rl_near))) / 12.0
        note = f"Supp. Fig. 6B filament count at [Nucleator]={conc} nM (rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(tv, row_c, y_axis, pc_x, half_r, x_axis=x_axis, slope_px=slope, note=note, clamp_y=(0.0, 1550.0))
        )
    return readings
