"""Calibrated reader for Fig. 1(b): 2-peak steady-state velocity profile v/U(x/l)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1 import read_curve


def read(image, at):
    return read_curve(image, at, "v")
