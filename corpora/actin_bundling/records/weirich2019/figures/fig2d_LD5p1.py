"""Calibrated reader for Fig. 2D — initial centering speed for L_D = 5.1 um (green diamonds)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2d import read_series


def read(image, at):
    return read_series(image, "LD5p1")
