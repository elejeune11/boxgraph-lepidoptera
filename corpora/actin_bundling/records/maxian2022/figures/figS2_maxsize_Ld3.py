"""Calibrated figure reader for figS2_maxsize_Ld3 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 60.0), y_bounds=(1.0, 1000.0), rgb=C.MATLAB_ORANGE, y_scale="log10", x_interior=[10,20,30,40,50], y_interior=[10.0,100.0], max_ink_gap=80)
