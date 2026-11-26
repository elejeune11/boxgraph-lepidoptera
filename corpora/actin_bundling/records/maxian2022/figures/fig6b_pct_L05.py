"""Calibrated figure reader for fig6b_pct_L05 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 30.0), y_bounds=(0.0, 100.0), rgb=C.MATLAB_PURPLE, x_interior=[10,20], y_interior=[20,40,60,80], exclude_fracs=[(0.20, 0.72, 0.98, 0.98)])
