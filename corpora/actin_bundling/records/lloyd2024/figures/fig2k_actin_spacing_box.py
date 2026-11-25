"""Digitize Fig. 2K: 50% developing H. sara F-actin bundle spacing notched boxplots (black x=1, blue x=2), p. 4, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_fig2_boxplot


def read(image, at):
    return read_fig2_boxplot(image, [1.6, 1.4, 1.2, 1.0, 0.8], "actin bundle spacing")
