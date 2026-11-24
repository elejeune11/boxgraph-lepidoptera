"""Calibrated reader for Fig. 4d (fig4d_alpha144)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2_4 import read_panel_curve

def read(image, at):
    return read_panel_curve(image, at, '4d', 'blue', style='solid', branch='top', dashed_x_range=(0.0, 10.6))
