"""Calibrated figure reader for fig7_bd_nb_two (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 5.0), y_bounds=(0.0, 5.0), rgb=C.FIG7_PINK, x_interior=[1,2,3,4], y_interior=[1,2,3,4], exclude_fracs=[(0.15, 0.01, 0.99, 0.22)], black_mode=False, tol=45.0, stem_search_px=6)
