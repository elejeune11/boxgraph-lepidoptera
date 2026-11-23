"""Calibrated reader for Main Fig. 3B left (p. 4): Src42 mutant ring spacing probability distribution."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3b_dist_wt import calibrate_fig3b_dist, read_circle_series


def read(image, at):
    x_axis, y_axis, box = calibrate_fig3b_dist(image)
    rgb = image.astype(float)
    green = (rgb[:, :, 1] > 180) & (rgb[:, :, 0] < 80) & (rgb[:, :, 2] < 80)
    return read_circle_series(green, x_axis, y_axis, box, at)
