"""Calibrated reader for Fig. 5(d): normalized mean tension vs density for S_0 = 0.1 (blue markers + 2-SD bars)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5 import read_fig5d_series


def read(image, at):
    return read_fig5d_series(image, "S0_0p1")
