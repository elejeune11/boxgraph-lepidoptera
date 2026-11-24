"""Calibrated reader for Fig. 10(a) unequal-core orange nullcline (dR_2/dt = 0) R_2/w vs R_1/w."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig10


def read(image, at):
    return _fig10.read_fig10a_nullcline(image, at, "R2")
