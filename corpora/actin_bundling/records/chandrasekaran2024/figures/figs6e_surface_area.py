"""Calibrated bar + SD error bar reader for Supplementary Fig. 6E actin-covered surface area fraction (1200 dpi)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

NUC_CONCS = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    y_axis = pc.fit_axis(info["y_ticks"], [1.0, 0.5, 0.0], max_residual_px=4.0)
    if len(info["x_ticks"]) != 8:
        raise ValueError(f"Expected 8 x-ticks in Supp. Fig. 6E, got {len(info['x_ticks'])}")
    c0, r0, c1, r1 = info["box"]

    readings = []
    for nuc, xt in zip(NUC_CONCS, info["x_ticks"]):
        cx = int(round(xt))
        # Left and right shoulders of green bar away from central red SD error bar
        bar_strip = np.concatenate(
            [image[r0 + 5:r1 - 3, cx - 60:cx - 22], image[r0 + 5:r1 - 3, cx + 22:cx + 60]],
            axis=1,
        ).astype(np.float64)
        # Green fill is ~RGB(103, 141, 65): use 50% transition of (G - B) ~ 38 for unbiased bar top
        greenness = np.median(bar_strip[:, :, 1] - bar_strip[:, :, 2], axis=1)
        g_rows = np.flatnonzero(greenness > 38.0) + (r0 + 5)
        if not g_rows.size:
            raise ValueError(f"Could not find green bar at [Nucleator]={nuc:g} nM in Supp. Fig. 6E")
        row_top = float(g_rows.min())

        eb_strip = image[r0 + 5:r1 - 3, cx - 16:cx + 16].astype(np.float64)
        redness = eb_strip[:, :, 0] - np.maximum(eb_strip[:, :, 1], eb_strip[:, :, 2])
        r_rows = np.flatnonzero(redness.max(axis=1) > 45.0) + (r0 + 5)
        eb_top = float(r_rows.min()) if r_rows.size else None
        eb_bot = float(r_rows.max() + 1) if r_rows.size else None

        # Panel is ~192 dpi rasterized (~6.25 px/src px at 1200 dpi): use honest +-4.5 px bar & +-5.5 px error-bar margins
        note = f"Supp. Fig. 6E surface area fraction at [Nucleator]={int(nuc)} nM (ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(
                nuc, row_top, y_axis, xt, half_px=4.5, note=note,
                bar_top_px=eb_top, bar_bot_px=eb_bot, bar_err_px=5.5, clamp_y=(0.0, 1.0),
            )
        )
    return readings
