"""Digitize SI Fig. 12C: 6 individual control ridge height dots (ordered ascending), p. 36, 1200 dpi."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs12 import read_s12c_dots


def read(image, at):
    return read_s12c_dots(image, group_index=1)
