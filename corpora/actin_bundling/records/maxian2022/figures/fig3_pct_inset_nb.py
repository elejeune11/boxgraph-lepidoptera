"""Calibrated figure reader for fig3_pct_inset_nb (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 12.0), y_bounds=(0.0, 55.0), rgb=C.MATLAB_BLUE, x_interior=[5,10], y_interior=[25,50])
