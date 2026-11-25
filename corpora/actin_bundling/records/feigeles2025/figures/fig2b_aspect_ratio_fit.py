"""Calibrated reader for Fig. 2B red exponential fit curve (tau ~ 5.8 s)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2b import read_fit


def read(image, at):
    return read_fit(image, at)
