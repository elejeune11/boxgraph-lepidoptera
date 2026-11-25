"""Calibrated reader for Fig. S5 — transverse free energy for L_M/L_D = 0.1 (green curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS5 import read_curve


def read(image, at):
    return read_curve(image, at, "ratio0p1")
