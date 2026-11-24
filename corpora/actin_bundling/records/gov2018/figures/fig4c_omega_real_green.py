"""Calibrated reader for Fig. 4(c) top: omega_real(q) small A* (solid green)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4c import read_fig4c


def read(image, at):
    return read_fig4c(image, at, 'omega_real_green')
