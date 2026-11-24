"""Calibrated reader for Supplementary Fig. 13C droplet-tip actin concentration vs aspect ratio (1200 dpi)."""

from __future__ import annotations
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import _plot_common as pc


def read(image: np.ndarray, at: list[float]) -> list[dict]:
    info = pc.find_plot_box(image)
    x_axis = pc.fit_axis(info["x_ticks"], [1.0, 2.0, 3.0, 4.0, 5.0, 5.56])
    y_axis = pc.fit_axis(info["y_ticks"], [80.0, 60.0, 40.0, 20.0, 0.0])
    c0, r0, c1, r1 = info["box"]

    # Pure black curve (grey reference lines have luminance ~190-220, so gray < 70 isolates black curve + text labels)
    gray = image.astype(np.float64).mean(axis=2)
    mask = np.zeros(gray.shape, dtype=bool)
    mask[r0 + 8:r1 - 8, c0 + 8:c1 - 8] = gray[r0 + 8:r1 - 8, c0 + 8:c1 - 8] < 70.0

    # Mask out the three text annotations ('C_total', 'C_total/2', 'C_total/3') which sit strictly above the curve
    # Locate their approximate pixel regions relative to calibrated coordinates so we never catch text strokes
    c_txt_1_lo, c_txt_1_hi = int(x_axis.to_pixel(2.2)), int(x_axis.to_pixel(3.4))
    r_txt_1_lo, r_txt_1_hi = int(y_axis.to_pixel(42.0)), int(y_axis.to_pixel(15.0))
    mask[r_txt_1_lo:r_txt_1_hi, c_txt_1_lo:c_txt_1_hi] = False

    c_txt_2_lo, c_txt_2_hi = int(x_axis.to_pixel(4.0)), int(x_axis.to_pixel(5.4))
    r_txt_2_lo, r_txt_2_hi = int(y_axis.to_pixel(16.0)), int(y_axis.to_pixel(2.5))
    mask[r_txt_2_lo:r_txt_2_hi, c_txt_2_lo:c_txt_2_hi] = False

    return pc.sample_curve_columns(
        mask, x_axis, y_axis, at, search_half_w=3, prefer="bottom",
        clamp_y=(0.0, 95.0), note_prefix="Supp. Fig. 13C tip actin concentration",
    )
