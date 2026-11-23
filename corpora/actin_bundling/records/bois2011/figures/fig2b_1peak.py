"""Calibrated reader for Fig. 2(b): Pe = 25, 1-peak steady-state potential (dashed)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2b import read_curve


def read(image, at):
    return read_curve(image, at, "1peak")
