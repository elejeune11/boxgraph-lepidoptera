"""Digitize Stage 16 (teal) aspect-ratio kernel density curve in Fig. 1f (PDF p. 2, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _density_common import read_density_curve


def read(image, at):
    return read_density_curve(
        image, at,
        x_labels=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0],
        y_labels_top_to_bottom=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0],
        curve="st16_teal",
    )
