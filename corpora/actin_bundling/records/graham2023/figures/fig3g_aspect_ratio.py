"""Calibrated reader for Main Fig. 3g (Shell, Ring, Disc aspect ratio)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _boxwhisker import read_boxwhisker


def read(image, at):
    return read_boxwhisker(image, [2.0,1.8,1.6,1.4,1.2,1.0], [1.0,2.0,3.0], y_max_search=1.70, label="Fig3g AR")
