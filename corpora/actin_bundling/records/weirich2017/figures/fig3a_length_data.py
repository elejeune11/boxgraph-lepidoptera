"""Calibrated reader for Fig. 3A measured average tactoid length vs time."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_length_data


def read(image, at):
    return read_length_data(image)
