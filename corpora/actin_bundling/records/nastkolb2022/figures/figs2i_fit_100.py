"""Digitize Fig. S2I log-normal fit curve at 1 uM VASP (purple curve; PDF p. 16, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs2i_common import read_figs2i_fit_100

def read(image, at):
    return read_figs2i_fit_100(image, at)
