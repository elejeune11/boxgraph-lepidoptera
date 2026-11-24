"""Calibrated reader for Fig. 9 rho=0.3 coupling path (neg)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig8_10 import read_fig9_path_rho03

def read(image, at):
    return read_fig9_path_rho03(image, at, 'neg')
