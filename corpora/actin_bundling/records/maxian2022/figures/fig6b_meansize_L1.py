"""Calibrated figure reader for fig6b_meansize_L1 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 30.0), y_bounds=(0.0, 40.0), rgb=C.MATLAB_YELLOW, x_interior=[10,20], y_interior=[10,20,30])
