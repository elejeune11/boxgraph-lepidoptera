"""Calibrated reader for Fig. 4-S1(d) col 2: total perpendicular tension sigma_perp (Reduce gamma by 200-fold)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4s1 import read_col2_tension


def read(image, at):
    return read_col2_tension(image, at, 'd', "perp")
