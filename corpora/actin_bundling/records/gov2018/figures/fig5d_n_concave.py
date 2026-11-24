"""Calibrated reader for Fig. 5(d) bottom: concave CMC coverage n_+(x)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5d import read_fig5d


def read(image, at):
    return read_fig5d(image, at, 'n_concave')
