"""Calibrated reader for Fig. 3b (fig3b_alpha16)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2_4 import read_panel_curve

def read(image, at):
    return read_panel_curve(image, at, '3b', 'blue', style='solid', branch='top', dashed_x_range=None)
