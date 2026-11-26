"""Calibrated figure reader for fig6a_bd_L05 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.0, 10.0), rgb=C.MATLAB_ORANGE, x_interior=[20,40], y_interior=[2,4,6,8])
