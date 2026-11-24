"""Calibrated reader for Fig. 4(B) right: High actin simulation curvature histogram (15 bins)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4 import read_hist_bars

def read(image, at):
    return read_hist_bars(image, 15)
