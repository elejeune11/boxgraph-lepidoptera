"""Calibrated figure reader for fig4_lm02_kon1 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 30.0), y_bounds=(0.0, 6.0), rgb=C.MATLAB_BLUE, x_interior=[5,10,15,20,25], y_interior=[1,2,3,4,5], exclude_fracs=[(0.58, 0.02, 0.99, 0.34)])
