"""Calibrated reader for Fig. 2C Inset — exponential fit curve d(t)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2c_inset import read_fit


def read(image, at):
    return read_fit(image, at)
