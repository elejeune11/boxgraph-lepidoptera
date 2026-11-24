"""Calibrated reader for SI Fig. S2D large-omega scaling curve (2*omega^(1/2))."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS2d import read_series


def read(image, at):
    return read_series(image, at, "large_w")
