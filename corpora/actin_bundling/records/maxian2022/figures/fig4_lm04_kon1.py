"""Calibrated figure reader for fig4_lm04_kon1 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 30.0), y_bounds=(0.0, 1.5), rgb=C.MATLAB_BLUE, x_interior=[5,10,15,20,25], y_interior=[0.5,1.0], max_ink_gap=80)
