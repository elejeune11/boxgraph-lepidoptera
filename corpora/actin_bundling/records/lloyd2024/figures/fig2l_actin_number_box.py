"""Digitize Fig. 2L: 50% developing H. sara F-actin bundle number notched boxplots (black x=1, blue x=2), p. 4, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_fig2_boxplot


def read(image, at):
    return read_fig2_boxplot(image, [50.0, 40.0, 30.0], "actin bundle number", has_dashed_lines=True)
