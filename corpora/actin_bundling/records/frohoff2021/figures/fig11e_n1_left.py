"""Calibrated reader for Fig. 11e (fig11e_n1_left)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig11_13 import read_fig11_13_curve

def read(image, at):
    return read_fig11_13_curve(image, at, '11e', 'green', style='solid', y_min_cut=0.78)
