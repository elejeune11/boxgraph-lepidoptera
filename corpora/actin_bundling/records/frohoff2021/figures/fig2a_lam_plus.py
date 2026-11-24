"""Calibrated reader for Fig. 2a (fig2a_lam_plus)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2_4 import read_panel_curve

def read(image, at):
    return read_panel_curve(image, at, '2a', 'blue', style='solid', branch='top', dashed_x_range=None)
