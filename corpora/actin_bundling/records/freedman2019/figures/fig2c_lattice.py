"""Calibrated reader for Fig. 2C — 1D equilibrium lattice model domain length vs Delta l_xl."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _panels


def read(image, at):
    return _panels.read_fig2c(image, at, "lattice")
