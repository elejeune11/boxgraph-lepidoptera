"""Calibrated reader for Fig. 3C red exponential fit curve (tau = 14.4 s)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3c import read_fit


def read(image, at):
    return read_fit(image, at)
