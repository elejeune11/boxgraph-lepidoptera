"""Digitize Fig. 2D log-normal fit curve at 0.125 uM VASP (red curve; PDF p. 4, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig2d_common import read_fig2d_fit

def read(image, at):
    return read_fig2d_fit(image, "0125", at)
