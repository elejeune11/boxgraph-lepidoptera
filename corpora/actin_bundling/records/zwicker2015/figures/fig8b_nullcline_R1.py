"""Calibrated reader for Fig. 8(b) autocatalytic two-droplet blue nullcline (Q_i = 300 w^3 k_0) R_2/w vs R_1/w."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig8


def read(image, at):
    return _fig8.read_fig8_nullcline(image, at, "R1")
