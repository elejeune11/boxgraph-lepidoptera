"""Calibrated reader for Fig. 1(d) solid black linear stability boundary P_c(beta)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1d import read_linear_boundary


def read(image, at):
    return read_linear_boundary(image, at)
