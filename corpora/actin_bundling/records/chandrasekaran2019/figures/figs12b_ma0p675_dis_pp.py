"""Calibrated reader for S12(b) (right, alpha:A=0.05, M:A=0.675, Dis++)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs12b

def read(image, at):
    return _figs12b.read_s12b_curve(image, 2, 'M:A=0.675 Dis++', at)
