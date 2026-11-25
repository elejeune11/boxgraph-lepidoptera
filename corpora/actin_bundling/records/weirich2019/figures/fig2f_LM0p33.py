"""Calibrated reader for Fig. 2F — continuum centering free energy for L_M = 0.33 um (red curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2f import read_curve


def read(image, at):
    return read_curve(image, at, "LM0p33")
