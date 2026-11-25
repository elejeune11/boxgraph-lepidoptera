"""Calibrated reader for Fig. 4-S1(a) col 2: total perpendicular tension sigma_perp (kappa = -0.2)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4s1 import read_col2_tension


def read(image, at):
    return read_col2_tension(image, at, 'a', "perp")
