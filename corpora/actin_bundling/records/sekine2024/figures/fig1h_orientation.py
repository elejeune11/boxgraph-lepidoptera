"""Digitize Actin orientation (Dominant direction vs Time, mean +/- SD) in Fig. 1h (PDF p. 2, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _timecourse_common import read_timecourse_ribbon


def read(image, at):
    return read_timecourse_ribbon(
        image, at,
        x_labels=[-80.0, -60.0, -40.0, -20.0, 0.0, 20.0],
        y_labels_top_to_bottom=[90.0, 60.0, 30.0, 0.0],
        color_kind="green",
    )
