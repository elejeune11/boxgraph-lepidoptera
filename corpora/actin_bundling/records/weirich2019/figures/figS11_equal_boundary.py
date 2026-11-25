"""Calibrated reader for Fig. S11 — elongated equal-division stability boundary (black dots)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS11 import read_boundary


def read(image, at):
    return read_boundary(image, at, "equal")
