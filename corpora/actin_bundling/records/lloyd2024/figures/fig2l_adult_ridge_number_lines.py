"""Digitize Fig. 2L: horizontal dashed adult cover ridge number reference lines (black x=1, blue x=2), p. 4, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_fig2l_adult_lines


def read(image, at):
    return read_fig2l_adult_lines(image)
