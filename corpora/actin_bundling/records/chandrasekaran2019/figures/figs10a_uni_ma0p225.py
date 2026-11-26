"""Calibrated reader for S10(a) (unipolar, alpha:A=0.1, M:A=0.225)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs10

def read(image, at):
    return _figs10.read_s10_series(image, 'unipolar', 0.225, at, read_sd=False)
