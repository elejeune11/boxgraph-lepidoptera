"""Calibrated reader for Fig. 2D — 1D equilibrium lattice model domain length vs L_p."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig2d(image, at, "lattice")
