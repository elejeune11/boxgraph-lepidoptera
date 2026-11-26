"""Calibrated figure reader for figS5_fit_b (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_black_fit(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.0, 2.0), x_interior=[10,20,30,40,50], y_interior=[0.5,1.0,1.5], branch="only", exclude_fracs=[(0.58, 0.67, 0.99, 0.99)])
