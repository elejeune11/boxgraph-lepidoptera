"""Calibrated reader for Fig. 2b (fig2b_rho00)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2_4 import read_panel_curve

def read(image, at):
    return read_panel_curve(image, at, '2b', 'green', style='solid', branch='top', dashed_x_range=None)
