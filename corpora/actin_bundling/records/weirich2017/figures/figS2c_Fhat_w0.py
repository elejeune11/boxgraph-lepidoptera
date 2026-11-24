"""Calibrated reader for SI Fig. S2C (omega = 0 blue curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS2c import read_curve


def read(image, at):
    return read_curve(image, at, "w0")
