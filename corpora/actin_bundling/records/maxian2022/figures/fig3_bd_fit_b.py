"""Calibrated figure reader for fig3_bd_fit_b (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_black_fit(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.0, 6.0), x_interior=[10,20,30,40,50], y_interior=[1,2,3,4,5], branch="lower", branch_split_x=15.5, exclude_fracs=[(0.57, 0.01, 0.995, 0.39)])
