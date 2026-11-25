"""Calibrated reader for Fig. S6 — harmonic spring constant k/F_0 for L_M = 1.5 um (blue curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6 import read_curve


def read(image, at):
    return read_curve(image, at, "LM1p5")
