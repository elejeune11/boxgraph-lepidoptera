"""Digitize Fig. S2I bundle width histogram markers at 1 uM VASP (purple diamonds; PDF p. 16, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs2i_common import read_figs2i_hist_100

def read(image, at):
    return read_figs2i_hist_100(image, at)
