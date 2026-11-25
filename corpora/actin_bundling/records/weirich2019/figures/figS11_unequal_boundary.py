"""Calibrated reader for Fig. S11 — elongated unequal-to-single transition boundary (red dots)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS11 import read_boundary


def read(image, at):
    return read_boundary(image, at, "unequal")
