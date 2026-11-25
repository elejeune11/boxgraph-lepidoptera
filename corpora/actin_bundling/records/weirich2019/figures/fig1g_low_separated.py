"""Calibrated reader for Fig. 1G Lower — separated morphology probability vs myosin concentration."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1g import read_series


def read(image, at):
    return read_series(image, "lower", "separated")
