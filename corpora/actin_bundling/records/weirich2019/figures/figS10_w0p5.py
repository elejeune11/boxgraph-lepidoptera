"""Calibrated reader for Fig. S10 — elongated two-tactoid free energy for wa^2*gamma/K^2 = 0.5 (red dashed curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS10 import read_curve


def read(image, at):
    return read_curve(image, at, "w0p5")
