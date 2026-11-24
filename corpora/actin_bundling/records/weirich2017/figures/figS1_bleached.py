"""Calibrated reader for SI Fig. S1 bleached-side mean fluorescence intensity."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1 import read_series


def read(image, at):
    return read_series(image, at, "bleached")
