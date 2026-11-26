"""Calibrated figure reader for fig2_link_k0007 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.0, 60.0), rgb=C.MATLAB_ORANGE, x_interior=[10,20,30,40,50], y_interior=[10,20,30,40,50], exclude_fracs=[(0.52, 0.67, 0.99, 0.98)], max_ink_gap=70)
