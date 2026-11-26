"""Calibrated figure reader for fig7_pct_b_one (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_standard_errorbar(image, at, x_bounds=(0.0, 5.0), y_bounds=(0.0, 60.0), rgb=C.FIG7_DGREEN, x_interior=[1,2,3,4], y_interior=[10,20,30,40,50], black_mode=False, tol=45.0, stem_search_px=6)
