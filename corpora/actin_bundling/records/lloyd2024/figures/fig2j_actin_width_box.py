"""Digitize Fig. 2J: 50% developing H. sara F-actin bundle width notched boxplots (black x=1, blue x=2), p. 4, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_fig2_boxplot


def read(image, at):
    return read_fig2_boxplot(image, [0.500, 0.475, 0.450, 0.425, 0.400], "actin bundle width")
