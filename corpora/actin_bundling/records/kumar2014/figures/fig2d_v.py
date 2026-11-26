"""Calibrated reader for Fig. 2(d) dimensionless flow velocity profile v*tau/ell at time t+T/4."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2 import read_fig2_profile


def read(image, at):
    return read_fig2_profile(image, at, "v", "d")
