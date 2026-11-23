"""Digitize Tube expansion (Tube diameter vs Time, mean and range) in Fig. 1j (PDF p. 2, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _timecourse_common import read_timecourse_ribbon


def read(image, at):
    return read_timecourse_ribbon(
        image, at,
        x_labels=[-80.0, -60.0, -40.0, -20.0, 0.0, 20.0],
        y_labels_top_to_bottom=[4.5, 4.0, 3.5, 3.0, 2.5],
        color_kind="blue",
    )
