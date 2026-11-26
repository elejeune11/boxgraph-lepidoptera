"""Calibrated figure reader for figS7_lm08_L1 (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 2.0), y_bounds=(0.2, 1.0), rgb=C.FIGS7_BLACK, x_interior=[0.5,1.0,1.5], y_interior=[0.3,0.4,0.5,0.6,0.7,0.8,0.9], exclude_fracs=[(0.17, 0.02, 0.99, 0.23)], black_mode=True, tol=45.0)
