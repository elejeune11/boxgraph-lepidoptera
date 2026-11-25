"""Calibrated reader for Fig. 3C measured coalescing condensate length vs time (black circles)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3c import read_data


def read(image, at):
    return read_data(image)
