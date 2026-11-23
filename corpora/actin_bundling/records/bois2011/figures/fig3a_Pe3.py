"""Calibrated reader for Fig. 3(a): Pe = 3 dispersion curve (bottom)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_curve


def read(image, at):
    return read_curve(image, at, 3, "Pe=3")
