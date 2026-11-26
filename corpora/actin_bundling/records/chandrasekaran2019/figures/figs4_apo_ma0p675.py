"""Calibrated reader for S4 Figure (apolar, M:A=0.675)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs4

def read(image, at):
    return _figs4.read_series(image, 'apolar', 0.675, at)
