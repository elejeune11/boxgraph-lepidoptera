"""Calibrated reader for Fig. 3a (fig3a_lam_plus)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2_4 import read_panel_curve

def read(image, at):
    return read_panel_curve(image, at, '3a', 'blue', style='solid', branch='top', dashed_x_range=None)
