"""Calibrated reader for Fig. 3C (p. 3): measured axial SAXS markers at R* = 0.1 (blue open squares)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3c import read_fig3c_data


def read(image, at):
    return read_fig3c_data(image, at, "R0p1")
