"""Calibrated reader for Fig. 4(c) bottom: omega_imag(q) small A* onset arch (solid green, left axis)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4c import read_fig4c


def read(image, at):
    return read_fig4c(image, at, 'omega_imag_green')
