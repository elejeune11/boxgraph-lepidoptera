"""Calibrated reader for Fig. 4-S1(b) col 3: active perpendicular tension (kappa = -0.5)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4s1 import read_col3_active


def read(image, at):
    return read_col3_active(image, at, 'b', "perp")
