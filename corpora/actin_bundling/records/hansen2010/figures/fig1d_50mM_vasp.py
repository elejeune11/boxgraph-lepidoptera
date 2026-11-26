"""Calibrated reader for Fig. 1D (p. 3): 50 mM KCl, 0-100 nM VASP bars with SD error bars."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig1d_series(image, '50mM_vasp')
