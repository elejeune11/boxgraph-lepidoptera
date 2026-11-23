"""Digitize mean longitudinal (0-30 deg, light-pink box black diamond) fusion duration vs time in Fig. 2b (PDF p. 4, 900 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fusion_boxplot_common import read_fusion_means


def read(image, at):
    return read_fusion_means(image, [-40.0, -30.0, -20.0, -10.0, 0.0], angle_slot=2)
