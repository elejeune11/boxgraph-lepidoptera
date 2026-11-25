"""Digitize Fig. 8N: all 20 individual cyto-D scale ridge-number dots, p. 11, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8 import read_8no_dots


def read(image, at):
    return read_8no_dots(image, y_labels=[60.0, 50.0, 40.0, 30.0], group_index=2, what="ridge number")
