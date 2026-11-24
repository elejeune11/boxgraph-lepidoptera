"""Calibrated reader for SI Fig. S2D small-omega scaling curve (1 + omega)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS2d import read_series


def read(image, at):
    return read_series(image, at, "small_w")
