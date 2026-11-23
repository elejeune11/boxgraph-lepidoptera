"""Digitize Stage 16 (teal) major-axis kernel density curve in Supplementary Fig. 1b (PDF p. 15, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _density_common import read_density_curve


def read(image, at):
    return read_density_curve(
        image, at,
        x_labels=[0.0, 0.5, 1.0, 1.5, 2.0],
        y_labels_top_to_bottom=[5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
        curve="st16_teal",
    )
