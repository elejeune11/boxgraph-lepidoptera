"""Calibrated figure reader for fig4_lm08_inset_kon4 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 15.0), y_bounds=(0.2, 0.45), rgb=C.MATLAB_YELLOW, x_interior=[5,10], y_interior=[0.3,0.4])
