"""Calibrated reader for Fig. 11a (fig11a_n1_minus)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig11_13 import read_fig11_13_curve

def read(image, at):
    return read_fig11_13_curve(image, at, '11a', 'green', style='dashed', y_min_cut=0.04)
