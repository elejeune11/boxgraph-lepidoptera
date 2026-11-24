"""Shared calibrated boxplot reader (orange median + 95% CI black whiskers) for Supplementary Figs. 5 and 10."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc


def read_ellipsoid_boxplot(image: np.ndarray, y_tick_vals: list[float], cat_labels: list[str],
                           panel_label: str) -> list[dict]:
    """Read orange medians and black 95% CI whisker ends from a boxplot panel in Supp. Fig. 5 or 10."""
    info = pc.find_plot_box(image)
    c0, r0, c1, r1 = info["box"]
    x_ticks = [t for t in info["x_ticks"] if c0 + 25 < t < c1 - 25]
    if len(x_ticks) != len(cat_labels):
        raise ValueError(f"{panel_label}: expected {len(cat_labels)} x-ticks, got {len(x_ticks)}")
    y_axis = pc.fit_axis(info["y_ticks"], y_tick_vals, max_residual_px=5.0)

    readings = []
    for idx, (lbl, xt) in enumerate(zip(cat_labels, x_ticks), start=1):
        cx = int(round(xt))
        strip = image[r0 + 5:r1 - 5, max(c0, cx - 80):min(c1, cx + 80)].astype(int)
        # Orange median line: R > 210, 95 < G < 175, B < 55
        is_org = (strip[:, :, 0] > 210) & (strip[:, :, 1] > 95) & (strip[:, :, 1] < 175) & (strip[:, :, 2] < 55)
        orows = np.flatnonzero(is_org.sum(axis=1) > 25) + (r0 + 5)
        if not orows.size:
            raise ValueError(f"{panel_label}: could not find orange median for {lbl}")
        r_med = 0.5 * (float(orows.min()) + float(orows.max())) + 0.5
        half_med = 0.5 * (float(orows.max()) - float(orows.min()) + 1.0)

        # Colored box interior (to locate box top/bottom edges for whisker stem walk)
        colored = ((strip.max(axis=2) - strip.min(axis=2)) > 25) & (strip.max(axis=2) > 80)
        brows = np.flatnonzero(colored.sum(axis=1) > 35) + (r0 + 5)
        box_top = int(brows.min()) if brows.size else int(orows.min())
        box_bot = int(brows.max()) if brows.size else int(orows.max())

        # Trace vertical black whisker stem upward from box_top and downward from box_bot near cx
        dark_col = (image[r0 + 5:r1 - 5, cx - 8:cx + 9].mean(axis=2) < 90).any(axis=1)
        H = len(dark_col)
        ru = max(0, box_top - (r0 + 5))
        while ru > 0 and (dark_col[ru - 1] or (ru > 1 and dark_col[ru - 2]) or (ru > 2 and dark_col[ru - 3])):
            ru -= 1
        rd = min(H - 1, box_bot - (r0 + 5))
        while rd < H - 1 and (dark_col[rd + 1] or (rd < H - 2 and dark_col[rd + 2]) or (rd < H - 3 and dark_col[rd + 3])):
            rd += 1

        w_top = float(ru + r0 + 5)
        w_bot = float(rd + r0 + 6)
        note = f"{panel_label} {lbl} median & 95% CI whiskers (ry={y_axis.residual_px:.2f}px)"
        readings.append(
            pc.make_reading(
                float(idx), r_med, y_axis, xt, half_px=half_med, note=note,
                bar_top_px=w_top, bar_bot_px=w_bot, bar_err_px=2.0, clamp_y=(0.0, 1.55),
            )
        )
    return readings
