"""Calibrated reader for SI Fig. S2 (left): Low actin simulation (dark blue circles)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS2 import read_colored_curve

def read(image, at):
    return read_colored_curve(image, at, "left", "low", branch="sim")
