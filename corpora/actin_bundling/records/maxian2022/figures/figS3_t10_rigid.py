"""Calibrated figure reader for figS3_t10_rigid (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS3_hist(image, at, x_bounds=(2.0, 17.0), y_bounds=(0.0, 25.0), rgb=C.MATLAB_YELLOW, x_interior=[5,10,15], y_interior=[5,10,15,20])
