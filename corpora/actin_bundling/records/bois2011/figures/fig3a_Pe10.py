"""Calibrated reader for Fig. 3(a): Pe = 10 dispersion curve (second from top)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_curve


def read(image, at):
    return read_curve(image, at, 1, "Pe=10")
