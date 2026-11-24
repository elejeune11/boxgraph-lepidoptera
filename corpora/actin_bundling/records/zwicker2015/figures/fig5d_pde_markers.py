"""Calibrated reader for Fig. 5(d) active PDE coarsening rate markers (orange discs) vs k_b/k_0."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5


def read(image, at):
    return _fig5.read_fig5d(image, at, "markers")
