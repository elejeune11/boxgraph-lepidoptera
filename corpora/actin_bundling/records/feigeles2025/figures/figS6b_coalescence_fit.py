"""Calibrated reader for Fig. S6B red exponential fit curve."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6b import read_fit


def read(image, at):
    return read_fit(image, at)
