"""Calibrated reader for Fig. 1(d) nonlinear oscillatory-to-stationary phase boundary P(beta)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1d import read_osc_to_stat_boundary


def read(image, at):
    return read_osc_to_stat_boundary(image, at)
