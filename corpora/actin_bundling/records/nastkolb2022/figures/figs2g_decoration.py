"""Digitize Fig. S2G bilayer actin decoration fraction (mean +/- SD) vs VASP concentration (PDF p. 16, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs2gh_s4h_common import read_vasp_marker_errorbars

def read(image, at):
    return read_vasp_marker_errorbars(image, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0], n_points=5)
