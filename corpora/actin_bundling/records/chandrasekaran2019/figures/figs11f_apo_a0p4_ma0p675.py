"""Calibrated reader for S11(f) (apolar, alpha:A=0.4, M:A=0.675)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs11

def read(image, at):
    return _figs11.read_s11_series(image, 'f', 0.675, at)
