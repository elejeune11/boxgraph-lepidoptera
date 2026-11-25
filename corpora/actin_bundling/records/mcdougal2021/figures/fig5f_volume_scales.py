"""Calibrated reader for Fig. 5F (p. 5): individual cover scale relative volumes (14 unobstructed pink crosses)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5cf import read_scales

def read(image, at):
    return read_scales(image, y_tick_vals=[100.0, 50.0], series="pink", has_legend=True)
