"""Calibrated reader for Fig. 4D — spherical unequal-to-single transition boundary (red dots)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4d import read_boundary


def read(image, at):
    return read_boundary(image, at, "unequal")
