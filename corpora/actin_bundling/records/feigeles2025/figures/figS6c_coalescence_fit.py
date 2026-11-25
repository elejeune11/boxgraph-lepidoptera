"""Calibrated reader for Fig. S6C red linear fit line."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS6c import read_fit


def read(image, at):
    return read_fit(image, at)
