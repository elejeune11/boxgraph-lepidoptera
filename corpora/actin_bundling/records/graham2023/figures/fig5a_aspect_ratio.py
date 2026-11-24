"""Calibrated reader for Main Fig. 5a (1:20, 1:10, 1:7.5 droplet aspect ratio)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _boxwhisker import read_boxwhisker


def read(image, at):
    return read_boxwhisker(image, [14.0,12.0,10.0,8.0,6.0,4.0,2.0,0.0], [0.05,0.1,0.1333], y_max_search=10.5, label="Fig5a AR")
