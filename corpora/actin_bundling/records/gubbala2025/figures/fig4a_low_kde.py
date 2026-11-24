"""Calibrated reader for Fig. 4(A) right: Low actin simulation KDE curve P(kappa)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4 import read_hist_kde

def read(image, at):
    return read_hist_kde(image, at)
