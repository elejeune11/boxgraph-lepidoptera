"""Calibrated reader for Fig. 3(d) finite-system stable stationary radius R_bar/w vs k_f/k_0 (solid blue curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3


def read(image, at):
    return _fig3.read_fig3bd(image, at, "Rbar")
