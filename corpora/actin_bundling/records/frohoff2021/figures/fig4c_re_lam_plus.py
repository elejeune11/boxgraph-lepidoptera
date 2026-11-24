"""Calibrated reader for Fig. 4c (fig4c_re_lam_plus)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2_4 import read_panel_curve

def read(image, at):
    return read_panel_curve(image, at, '4c', 'green', style='solid', branch='top', dashed_x_range=(3.7, 9.9))
