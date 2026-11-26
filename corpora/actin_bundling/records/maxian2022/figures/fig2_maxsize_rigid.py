"""Calibrated figure reader for fig2_maxsize_rigid (maxian2022)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig2_maxsize(image, at, C.MATLAB_YELLOW)
