"""Digitize dotted Gaussian fit curve in Supplementary Fig. 7c2 (3-M labyrinth, PDF p. 22, 1200 dpi)."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _s7c_common import read_s7c_fit


def read(image, at):
    return read_s7c_fit(image, at, [0.0, 2.0, 4.0, 6.0], [0.8, 0.6, 0.4, 0.2, 0.0])
