"""Calibrated reader for Fig. 4-S1(b) col 2: total parallel tension sigma_par (kappa = -0.5)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4s1 import read_col2_tension


def read(image, at):
    return read_col2_tension(image, at, 'b', "par")
