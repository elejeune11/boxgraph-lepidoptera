"""Calibrated reader for Main Fig. 6b FRAP average recovery curve at 10% PEG (blue)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig6b import read_frap


def read(image, at):
    return read_frap(image, at, "10peg")
