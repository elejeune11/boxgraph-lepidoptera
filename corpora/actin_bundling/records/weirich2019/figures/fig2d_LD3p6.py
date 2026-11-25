"""Calibrated reader for Fig. 2D — initial centering speed for L_D = 3.6 um (black x's)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2d import read_series


def read(image, at):
    return read_series(image, "LD3p6")
