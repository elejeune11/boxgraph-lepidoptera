"""Calibrated reader for S11(b) (unipolar, alpha:A=0.1, M:A=0.0225)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs11

def read(image, at):
    return _figs11.read_s11_series(image, 'b', 0.0225, at)
