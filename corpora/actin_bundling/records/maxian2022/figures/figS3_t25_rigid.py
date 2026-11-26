"""Calibrated figure reader for figS3_t25_rigid (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS3_hist(image, at, x_bounds=(2.0, 158.0), y_bounds=(0.0, 30.0), rgb=C.MATLAB_YELLOW, x_scale="log10", x_interior=[4,6,10,18,31,53,92], y_interior=[5,10,15,20,25])
