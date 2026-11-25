"""Digitize Fig. 8N: ridge number mean +- SD for control (x=1) and cyto-D (x=2), p. 11, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8 import read_8no_mean_sd


def read(image, at):
    return read_8no_mean_sd(image, y_labels=[60.0, 50.0, 40.0, 30.0], what="ridge number")
