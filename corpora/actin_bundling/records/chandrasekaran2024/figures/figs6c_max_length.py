"""Calibrated reader for Supplementary Fig. 6C upper filament length (at CDF = 0.96) vs [Nucleator] (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

NUC_CONCS_LR = [150.0, 125.0, 100.0, 75.0, 50.0, 25.0, 10.0, 0.0]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read upper filament length (at CDF = 0.96) across all 8 [Nucleator] concentrations."""
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], max_residual_px=4.0)
    y_axis = pc.fit_axis(info["y_ticks"], [1.0, 0.75, 0.5, 0.25, 0.0], max_residual_px=4.0)
    c0, r0, c1, r1 = info["box"]

    arr = image.astype(int)
    colored = ((arr.max(axis=2) - arr.min(axis=2)) > 65) & (arr.min(axis=2) < 120)
    colored[:r0 + 5, :] = False
    colored[r1 - 5:, :] = False
    colored[:, :c0 + 5] = False
    colored[:, c1 - 5:] = False

    row_96 = y_axis.to_pixel(0.96)
    ri = int(round(row_96))
    cs = np.flatnonzero(colored[ri - 2:ri + 3, :].any(axis=0))
    sp = np.flatnonzero(np.diff(cs) > 6)
    st = np.concatenate(([0], sp + 1))
    en = np.concatenate((sp, [cs.size - 1]))
    runs = [(0.5 * (float(cs[a]) + float(cs[b])) + 0.5, 0.5 * (float(cs[b]) - float(cs[a]) + 1.0))
            for a, b in zip(st, en) if (cs[b] - cs[a] + 1) >= 5]
    if len(runs) != 8:
        raise ValueError(f"Expected 8 CDF curves at CDF=0.96 in Supp. Fig. 6C, found {len(runs)}")

    px_per_um = abs(x_axis.to_pixel(6.0) - x_axis.to_pixel(0.0)) / 6.0
    x_misfit = max(x_axis.residual_px, 0.5)

    readings = []
    # Return sorted by ascending [Nucleator] (0..150 nM)
    for conc, (col_c, half_c) in reversed(list(zip(NUC_CONCS_LR, runs))):
        L_val = x_axis.to_value(col_c)
        err_um = (half_c + 1.0 + x_misfit) / px_per_um
        readings.append({
            "x": float(conc),
            "y": float(L_val),
            "lo": float(max(0.0, L_val - err_um)),
            "hi": float(min(6.5, L_val + err_um)),
            "px": [float(col_c), float(row_96)],
            "px_lo": float(row_96 - 10.0),
            "px_hi": float(row_96 + 10.0),
            "note": f"Supp. Fig. 6C filament length at CDF=0.96 for [Nucleator]={int(conc)} nM (rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)",
        })
    return readings
