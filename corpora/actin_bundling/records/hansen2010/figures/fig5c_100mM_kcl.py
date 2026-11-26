"""Calibrated reader for Fig. 5C (p. 8): 2 uM profilin-actin at 100 mM KCl (0, 25, 50 nM VASP)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig5c_bars(image)
