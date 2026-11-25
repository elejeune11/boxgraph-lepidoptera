"""Calibrated reader for Fig. 2B measured droplet aspect ratio vs time (blue filled circles)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2b import read_data


def read(image, at):
    return read_data(image)
