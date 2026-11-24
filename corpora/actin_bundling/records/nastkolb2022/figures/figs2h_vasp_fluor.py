"""Digitize Fig. S2H bilayer VASP fluorescence intensity [AU] vs VASP concentration (PDF p. 16, 1200 dpi)."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figs2gh_s4h_common import read_vasp_marker_errorbars

def read(image, at):
    return read_vasp_marker_errorbars(image, [0.0, 0.5, 1.0, 1.5, 2.0], [600.0, 500.0, 400.0, 300.0, 200.0, 100.0, 0.0], n_points=6)
