"""Digitize Fig. S4F bundle width boxplots (1=VASP, 2=VASP noCP SLB, 3=VASP noCP ML; PDF p. 18, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs4f_common import read_figs4f

def read(image, at):
    return read_figs4f(image, "width")
