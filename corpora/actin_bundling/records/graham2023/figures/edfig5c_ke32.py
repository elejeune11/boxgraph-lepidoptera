"""Calibrated reader for Extended Data Fig. 5c: f(e) = 32(e - 1) (red curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _edfig5c import read_curve


def read(image, at):
    return read_curve(image, at, "ke32")
