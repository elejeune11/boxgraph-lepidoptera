"""Calibrated reader for Fig. S6C measured coalescence Delta L vs tau (6 black circles)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6c import read_data


def read(image, at):
    return read_data(image)
