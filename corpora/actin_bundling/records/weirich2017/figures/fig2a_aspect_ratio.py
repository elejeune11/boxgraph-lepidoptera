"""Calibrated reader for Fig. 2A tactoid aspect ratio L/r vs filamin %."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2a import read_series


def read(image, at):
    return read_series(image, "aspect_ratio")
