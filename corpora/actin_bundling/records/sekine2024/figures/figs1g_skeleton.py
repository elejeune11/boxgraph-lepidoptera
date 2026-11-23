"""Digitize Actin elongation (Mean skeleton length vs Time, mean +/- SD) in Supplementary Fig. 1g (PDF p. 15, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _timecourse_common import read_timecourse_ribbon


def read(image, at):
    return read_timecourse_ribbon(
        image, at,
        x_labels=[-80.0, -60.0, -40.0, -20.0, 0.0, 20.0],
        y_labels_top_to_bottom=[0.75, 0.50, 0.25],
        color_kind="green",
    )
