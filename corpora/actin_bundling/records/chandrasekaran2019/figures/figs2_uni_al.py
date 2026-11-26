"""Calibrated reader for S2 Figure (unipolar, blue AL cluster)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs2

def read(image, at):
    return _figs2.read_s2_al(image, 'unipolar', at)
