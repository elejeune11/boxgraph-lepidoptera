"""Calibrated figure reader for figS3_t5_rigid (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS3_hist(image, at, x_bounds=(2.0, 8.0), y_bounds=(0.0, 20.0), rgb=C.MATLAB_YELLOW, x_interior=[3,4,5,6,7], y_interior=[5,10,15], exclude_fracs=[(0.62, 0.04, 0.98, 0.34)])
