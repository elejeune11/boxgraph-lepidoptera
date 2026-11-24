"""Calibrated reader for Fig. 8d (fig8d_n1_minus)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8_10 import read_fig8_10_curve

def read(image, at):
    return read_fig8_10_curve(image, at, '8d', 'blue', style='dashed', dashed_above_x=None)
