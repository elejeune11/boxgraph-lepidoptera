"""Digitize Fig. 3F pore size boxplots (broken right y-axis) across 5 protein conditions (PDF p. 5, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3f_common import read_fig3f

def read(image, at):
    return read_fig3f(image, "pore")
