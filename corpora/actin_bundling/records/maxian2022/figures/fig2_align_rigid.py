"""Calibrated figure reader for fig2_align_rigid (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.4, 1.0), rgb=C.MATLAB_YELLOW, x_interior=[10,20,30,40,50], y_interior=[0.5,0.6,0.7,0.8,0.9])
