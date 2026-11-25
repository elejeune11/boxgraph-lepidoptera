"""Calibrated reader for Fig. 5(d): theoretical dashed line lambda = 1.0."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5 import read_fig5d_theory_line


def read(image, at):
    return read_fig5d_theory_line(image, at)
