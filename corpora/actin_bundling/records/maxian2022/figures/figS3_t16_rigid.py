"""Calibrated figure reader for figS3_t16_rigid (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS3_hist(image, at, x_bounds=(3.0, 46.0), y_bounds=(0.0, 40.0), rgb=C.MATLAB_YELLOW, x_interior=[5,7,9,11,14,18,22,26,30,34,38,42], y_interior=[10,20,30])
