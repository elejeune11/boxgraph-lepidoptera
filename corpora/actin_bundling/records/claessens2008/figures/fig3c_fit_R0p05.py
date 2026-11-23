"""Calibrated reader for Fig. 3C (p. 3): continuous four-sphere fit curve at R* = 0.05 (red)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3c import read_fig3c_fit


def read(image, at):
    return read_fig3c_fit(image, at, "R0p05")
