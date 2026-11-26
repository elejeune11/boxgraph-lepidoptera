"""Calibrated figure reader for fig3_align_nb (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.7, 1.0), rgb=C.MATLAB_BLUE, x_interior=[10,20,30,40,50], y_interior=[0.75,0.8,0.85,0.9,0.95])
