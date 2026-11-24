"""Digitize Fig. S2F raw TIRFM actin intensity kinetics at 0.5 uM VASP (green curve; PDF p. 16, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs2f_common import read_figs2f

def read(image, at):
    return read_figs2f(image, "050", at)
