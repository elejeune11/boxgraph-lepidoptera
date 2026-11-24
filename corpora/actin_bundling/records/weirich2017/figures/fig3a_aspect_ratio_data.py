"""Calibrated reader for Fig. 3A measured tactoid aspect ratio L/r vs time."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_aspect_ratio_data


def read(image, at):
    return read_aspect_ratio_data(image)
