"""Calibrated figure reader for figS4_lm04_kon025 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 2.0), y_bounds=(0.2, 0.9), rgb=C.MATLAB_ORANGE, x_interior=[0.5,1.0,1.5], y_interior=[0.3,0.4,0.5,0.6,0.7,0.8])
