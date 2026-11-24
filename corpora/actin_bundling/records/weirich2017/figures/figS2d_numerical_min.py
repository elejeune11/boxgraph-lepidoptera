"""Calibrated reader for SI Fig. S2D numerical minimization dots."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS2d import read_series


def read(image, at):
    return read_series(image, at, "numerical_min")
