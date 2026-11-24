"""Calibrated reader for Fig. 1C exponential fit curve (tau_R = 880 s)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1c import read_fit


def read(image, at):
    return read_fit(image, at)
