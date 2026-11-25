"""Calibrated reader for Fig. 5D (p. 5): 6-timepoint mean cover scale width curve."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5cf import read_means

def read(image, at):
    return read_means(image, y_tick_vals=[50.0], at=at, series="black", has_legend=False)
