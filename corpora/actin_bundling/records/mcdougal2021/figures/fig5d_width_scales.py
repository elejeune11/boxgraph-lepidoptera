"""Calibrated reader for Fig. 5D (p. 5): individual measured cover scale widths (15 crosses)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5cf import read_scales

def read(image, at):
    return read_scales(image, y_tick_vals=[50.0], series="black", has_legend=False)
