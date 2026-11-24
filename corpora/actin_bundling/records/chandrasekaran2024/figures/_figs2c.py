"""Shared calibrated bar + SD error bar reader for Supplementary Fig. 2C (900 dpi, 2 columns x 3 rows = 36 conditions)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc

K_UNBIND_VALS = [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0]


def read_figs2c(image: np.ndarray, row_idx: int, k_bind_label: str) -> list[dict]:
    """Read one row (0=top, 1=mid, 2=bottom) of a 3-panel column in Supp. Fig. 2C."""
    boxes = pc.find_stacked_boxes(image, 3)
    b = boxes[row_idx]
    c0, r0, c1, r1 = b["box"]
    x_ticks = [t for t in b["x_ticks"] if c0 + 25 < t < c1 - 25]
    if len(x_ticks) != 6:
        raise ValueError(f"Expected 6 x-ticks in Supp. Fig. 2C ({k_bind_label}), got {len(x_ticks)}")
    y_axis = pc.fit_axis(b["y_ticks"], [1.0, 0.5, 0.0], max_residual_px=3.0)

    readings = []
    for ku, xt in zip(K_UNBIND_VALS, x_ticks):
        cx = int(round(xt))
        strip = image[r0 + 4:r1 - 2, max(c0, cx - 75):min(c1, cx + 75)].astype(int)
        is_green = (strip[:, :, 1] - strip[:, :, 2] > 25) & (strip[:, :, 1] > 90)
        g_rows = np.flatnonzero(is_green.mean(axis=1) > 0.30) + (r0 + 4)
        if not g_rows.size:
            raise ValueError(f"Could not find green bar at k_bind={k_bind_label}, k_unbind={ku:g}")
        row_top = float(g_rows.min())

        is_red = (strip[:, :, 0] > 140) & (strip[:, :, 1] < 80) & (strip[:, :, 2] < 80)
        r_rows = np.flatnonzero(is_red.any(axis=1)) + (r0 + 4)
        eb_top = float(r_rows.min()) if r_rows.size else None
        eb_bot = float(r_rows.max() + 1) if r_rows.size else None

        note = f"Supp. Fig. 2C final surface area fraction at k_bind={k_bind_label}, k_unbind={ku:g}/s (ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(
                ku, row_top, y_axis, xt, half_px=1.5, note=note,
                bar_top_px=eb_top, bar_bot_px=eb_bot, bar_err_px=2.0, clamp_y=(0.0, 1.0),
            )
        )
    return readings
