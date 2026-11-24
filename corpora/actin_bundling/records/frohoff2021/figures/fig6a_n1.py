"""Calibrated reader for Fig. 6a (fig6a_n1)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5_6 import read_fig56_curve

def read(image, at):
    return read_fig56_curve(image, at, '6a', 'blue', style='solid', branch_idx=0, dashed_above_x=None)
