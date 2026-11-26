"""Calibrated reader for S12(b) (left, alpha:A=0.05, M:A=0.0225, Dis--)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs12b

def read(image, at):
    return _figs12b.read_s12b_curve(image, 0, 'M:A=0.0225 Dis--', at)
