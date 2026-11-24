"""Shared calibrated stacked-bar reader for Supplementary Fig. 11D high-actin network probabilities (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

VASP_CONCS = [2.0, 0.8, 0.4]
N_FILS = [30.0, 75.0, 120.0, 165.0]


def read_figs11d(image: np.ndarray, panel_idx: int, category: str) -> list[dict]:
    """Read one [VASP-tet] subpanel (0=2.0, 1=0.8, 2=0.4 uM) and one cluster category
    ('ring_shell_sb' green bottom, 'ring_multiring' orange middle, 'locally_bundled' red top)."""
    boxes = pc.find_stacked_boxes(image, 3)
    b = boxes[panel_idx]
    y_axis = pc.fit_axis(b["y_ticks"], [1.0, 0.8, 0.6, 0.4, 0.2, 0.0], max_residual_px=3.5)
    if len(b["x_ticks"]) != 4:
        raise ValueError(f"Expected 4 x-ticks in Supp. Fig. 11D panel {panel_idx}, got {len(b['x_ticks'])}")
    c0, r0, c1, r1 = b["box"]
    row_0 = y_axis.to_pixel(0.0)
    row_1 = y_axis.to_pixel(1.0)
    px_per_unit = abs(row_0 - row_1)
    err_val = (2.5 + max(y_axis.residual_px, 0.5)) / px_per_unit
    conc = VASP_CONCS[panel_idx]

    readings = []
    for nf, xt in zip(N_FILS, b["x_ticks"]):
        cx = int(round(xt))
        col_rgb = image[r0 + 2:r1 - 1, cx - 35:cx + 35].astype(np.float64).mean(axis=1)
        is_g = (col_rgb[:, 1] - col_rgb[:, 0] > 15.0) & (col_rgb[:, 1] > 130.0) & (col_rgb[:, 2] < 110.0)
        is_o = (col_rgb[:, 0] > 165.0) & (col_rgb[:, 1] > 95.0) & (col_rgb[:, 1] < 155.0) & (col_rgb[:, 2] < 70.0)
        is_r = (col_rgb[:, 0] > 125.0) & (col_rgb[:, 1] < 75.0) & (col_rgb[:, 2] < 75.0)

        rg = np.flatnonzero(is_g) + (r0 + 2)
        ro = np.flatnonzero(is_o) + (r0 + 2)
        rr = np.flatnonzero(is_r) + (r0 + 2)

        # Stack order bottom -> top: green (Ring-Shell w. strong bundling), orange (Ring/Multi-ring), red (Locally-bundled)
        has_g, has_o, has_r = rg.size >= 8, ro.size >= 8, rr.size >= 8
        if has_g:
            row_g_top = 0.5 * (float(rg.min()) + float(ro.max()) + 1.0) if has_o else (
                0.5 * (float(rg.min()) + float(rr.max()) + 1.0) if has_r else row_1
            )
            pg = min(1.0, max(0.0, y_axis.to_value(row_g_top)))
        else:
            row_g_top = row_0
            pg = 0.0

        if has_r:
            row_r_bot = 0.5 * (float(rr.max()) + 1.0 + float(ro.min())) if has_o else (
                0.5 * (float(rr.max()) + 1.0 + float(rg.min())) if has_g else row_0
            )
            pr = min(1.0, max(0.0, 1.0 - y_axis.to_value(row_r_bot)))
        else:
            row_r_bot = row_1
            pr = 0.0

        po = min(1.0, max(0.0, 1.0 - pg - pr)) if has_o else 0.0

        if category == "ring_shell_sb":
            val = pg
            px_row = row_g_top if val > 0 else row_0
            note = f"Supp. Fig. 11D [VASP-tet]={conc} uM Ring-Shell-SB (green) at N={int(nf)} (ry={y_axis.residual_px:.2f}px)"
        elif category == "ring_multiring":
            val = po
            px_row = row_r_bot if val > 0 else row_0
            note = f"Supp. Fig. 11D [VASP-tet]={conc} uM Ring/Multi-ring (orange) at N={int(nf)} (ry={y_axis.residual_px:.2f}px)"
        elif category == "locally_bundled":
            val = pr
            px_row = row_r_bot if val > 0 else row_1
            note = f"Supp. Fig. 11D [VASP-tet]={conc} uM Locally-bundled (red) at N={int(nf)} (ry={y_axis.residual_px:.2f}px)"
        else:
            raise ValueError(f"Unknown category {category!r}")

        h_px = err_val * px_per_unit
        readings.append({
            "x": float(nf),
            "y": float(val),
            "lo": float(max(0.0, val - err_val)),
            "hi": float(min(1.0, val + err_val)),
            "px": [float(xt), float(px_row)],
            "px_lo": float(px_row - h_px),
            "px_hi": float(px_row + h_px),
            "note": note,
        })
    return readings
