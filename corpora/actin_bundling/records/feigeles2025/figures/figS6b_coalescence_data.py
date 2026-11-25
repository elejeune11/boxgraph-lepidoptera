"""Calibrated reader for Fig. S6B measured free-solution coalescence length vs time (12 black circles)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6b import read_data


def read(image, at):
    return read_data(image)
