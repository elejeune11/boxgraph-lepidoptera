"""Calibrated reader for Fig. 1C measured FRAP recovery points."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1c import read_data


def read(image, at):
    return read_data(image)
