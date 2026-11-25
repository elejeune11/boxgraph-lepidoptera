"""Calibrated reader for Fig. 5(c) (S_0 = 0.3): solid time-averaged parallel tension."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5 import read_fig5c


def read(image, at):
    return read_fig5c(image, at, 'S0_0p3', 'par')
