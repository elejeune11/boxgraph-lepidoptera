"""Calibrated reader for Fig. 10a (fig10a_n3)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8_10 import read_fig8_10_curve

def read(image, at):
    return read_fig8_10_curve(image, at, '10a', 'orange', style='solid', dashed_above_x=None)
