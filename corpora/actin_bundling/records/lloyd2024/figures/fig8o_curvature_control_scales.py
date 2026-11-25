"""Digitize Fig. 8O: all 12 individual control scale curvature dots, p. 11, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8 import read_8no_dots


def read(image, at):
    return read_8no_dots(image, y_labels=[0.15, 0.10, 0.05, 0.00], group_index=1, what="ridge curvature")
