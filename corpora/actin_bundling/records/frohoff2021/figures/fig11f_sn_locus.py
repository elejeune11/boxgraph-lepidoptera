"""Calibrated reader for Fig. 11(f) (fig11f_sn_locus)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig11_13 import read_fig11f_locus

def read(image, at):
    return read_fig11f_locus(image, at, 'sn')
