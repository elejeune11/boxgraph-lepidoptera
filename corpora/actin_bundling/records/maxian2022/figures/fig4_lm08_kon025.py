"""Calibrated figure reader for fig4_lm08_kon025 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.0, 0.5), rgb=C.MATLAB_ORANGE, x_interior=[10,20,30,40,50], y_interior=[0.1,0.2,0.3,0.4], exclude_fracs=[(0.42, 0.03, 0.99, 0.57)])
