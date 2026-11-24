"""Calibrated reader for Supplementary Fig. 2B icosphere triangle-side boxplot (900 dpi)."""

from __future__ import annotations
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], max_residual_px=2.0)
    y_vals = [1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    y_axis = pc.fit_axis(info["y_ticks"], y_vals, max_residual_px=3.5)
    c0, r0, c1, r1 = info["box"]

    readings = []
    for n, xt in enumerate(info["x_ticks"]):
        cx = int(round(xt))
        sub = image[r0 + 5:r1 - 5, max(c0, cx - 110):min(c1, cx + 110)].astype(int)
        # Median is drawn as a horizontal blue bar
        is_blue = (sub[:, :, 2] - sub[:, :, 0] > 35) & (sub[:, :, 2] > 100)
        b_rows = np.flatnonzero(is_blue.sum(axis=1) > 20) + (r0 + 5)
        if not b_rows.size:
            raise ValueError(f"Could not find blue median line at subdivision={n} in Supp. Fig. 2B")
        row_med = 0.5 * (float(b_rows.min()) + float(b_rows.max())) + 0.5
        half_med = 0.5 * (float(b_rows.max()) - float(b_rows.min()) + 1.0)

        # 95% confidence range is drawn in red
        is_red = (sub[:, :, 0] - sub[:, :, 2] > 50) & (sub[:, :, 0] > 150) & (sub[:, :, 1] < 100)
        r_rows = np.flatnonzero(is_red.sum(axis=1) > 10) + (r0 + 5)
        if not r_rows.size:
            raise ValueError(f"Could not find red 95% CI range at subdivision={n} in Supp. Fig. 2B")
        ci_top = float(r_rows.min())
        ci_bot = float(r_rows.max() + 1)

        note = f"Supp. Fig. 2B blue median & red 95% CI at subdivision={n} (rx={x_axis.residual_px:.2f}px, ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(
                float(n), row_med, y_axis, xt, half_px=half_med, x_axis=x_axis, note=note,
                bar_top_px=ci_top, bar_bot_px=ci_bot, bar_err_px=2.0, clamp_y=(0.0, 1.2),
            )
        )
    return readings
