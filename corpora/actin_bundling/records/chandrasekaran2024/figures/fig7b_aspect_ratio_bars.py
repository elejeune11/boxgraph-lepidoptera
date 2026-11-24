"""Calibrated reader for the 10 simulation mean aspect-ratio bars and SD error bars in Fig. 7B (1200 dpi)."""

from __future__ import annotations
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))

import math
import numpy as np
import _plot_common as pc

FRACS = [2 / 3, 5 / 6, 11 / 12, 1.0, 13 / 12, 7 / 6, 5 / 4, 4 / 3, 5 / 3, 2.0]


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    x_vals = [f * math.pi for f in FRACS]
    x_axis = pc.fit_axis(info["x_ticks"], x_vals)
    y_axis = pc.fit_axis(info["y_ticks"], [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    c0, r0, c1, r1 = info["box"]

    readings = []
    for f, xv, xt in zip(FRACS, x_vals, info["x_ticks"]):
        cx = int(round(xt))
        strip = image[r0 + 3:r1 - 3, cx - 38:cx + 38].astype(int)
        is_green = (
            (strip[:, :, 1] - strip[:, :, 0] > 15)
            & (strip[:, :, 1] - strip[:, :, 2] > 25)
            & (strip[:, :, 1] > 90)
        )
        g_rows = np.flatnonzero(is_green.mean(axis=1) > 0.25) + (r0 + 3)
        if not g_rows.size:
            raise ValueError(f"Could not locate green bar at L={f:.4f}*pi in Fig. 7B")
        row_top = float(g_rows.min())

        # Red SD error bar centered near bar center
        is_red = (strip[:, :, 0] > 130) & (strip[:, :, 1] < 70) & (strip[:, :, 2] < 70)
        r_rows = np.flatnonzero(is_red.any(axis=1)) + (r0 + 3)
        if not r_rows.size:
            raise ValueError(f"Could not locate red SD error bar at L={f:.4f}*pi in Fig. 7B")
        eb_top = float(r_rows.min())
        eb_bot = float(r_rows.max() + 1)

        note = f"Fig. 7B mean aspect ratio bar at L_max={f:.4f}*pi*R_drop (rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(
                xv, row_top, y_axis, xt, half_px=1.5, x_axis=x_axis, note=note,
                bar_top_px=eb_top, bar_bot_px=eb_bot, bar_err_px=2.0,
            )
        )
    return readings
