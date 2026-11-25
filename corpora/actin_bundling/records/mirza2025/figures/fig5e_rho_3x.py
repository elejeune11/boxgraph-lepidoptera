"""Calibrated reader for Fig. 5(e): deviatoric tension ratio vs S_0 at 3*rho (red)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5 import read_fig5e_series


def read(image, at):
    return read_fig5e_series(image, 'rho_3x')
