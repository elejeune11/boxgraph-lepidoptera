"""Calibrated reader for Fig. 2(a): Pe = 3 potential curve (dotted)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2a import read_curve


def read(image, at):
    return read_curve(image, at, "Pe3")
