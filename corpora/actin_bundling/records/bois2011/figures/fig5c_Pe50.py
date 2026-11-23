"""Calibrated reader for Fig. 5(c): ASDM Pe = 50 dispersion curve at 'x'."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5c import read_curve


def read(image, at):
    return read_curve(image, at, "Pe50")
