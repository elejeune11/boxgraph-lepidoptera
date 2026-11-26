"""Calibrated reader for S4 Figure (unipolar, M:A=0.225)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs4

def read(image, at):
    return _figs4.read_series(image, 'unipolar', 0.225, at)
