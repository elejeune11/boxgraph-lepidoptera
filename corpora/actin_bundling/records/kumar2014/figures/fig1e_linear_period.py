"""Calibrated reader for Fig. 1(e) linear stability period curve 2*pi/Im[tau*lambda_+(k_1)]."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1e import read_linear_period


def read(image, at):
    return read_linear_period(image, at)
