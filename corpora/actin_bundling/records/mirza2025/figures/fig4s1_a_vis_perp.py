"""Calibrated reader for Fig. 4-S1(a) col 3: viscous perpendicular tension (kappa = -0.2)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4s1 import read_col3_viscous


def read(image, at):
    return read_col3_viscous(image, at, 'a', "perp")
