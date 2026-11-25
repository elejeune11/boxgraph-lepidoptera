"""Calibrated reader for Fig. 5(e): dashed linear regression line (kappa = -1.6)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5 import read_fig5e_fit_line


def read(image, at):
    return read_fig5e_fit_line(image, at)
