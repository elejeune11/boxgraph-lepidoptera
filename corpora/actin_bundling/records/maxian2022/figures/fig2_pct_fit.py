"""Calibrated figure reader for fig2_pct_fit (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_black_fit(image, at, x_bounds=(0.0, 60.0), y_bounds=(0.0, 100.0), x_interior=[10,20,30,40,50], y_interior=[20,40,60,80], branch="only", exclude_fracs=[(0.50, 0.72, 0.99, 0.98)], half_win_px=20)
