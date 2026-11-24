"""Calibrated reader for Fig. 5c (fig5c_phi2)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5_6 import read_fig56_curve

def read(image, at):
    return read_fig56_curve(image, at, '5c', 'orange', style='solid', branch_idx=0, dashed_above_x=None)
