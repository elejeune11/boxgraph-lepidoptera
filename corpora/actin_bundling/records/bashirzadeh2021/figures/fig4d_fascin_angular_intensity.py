"""Calibrated reader for Fig. 4d: Fascin fluorescence intensity (a.u., left y-axis 1400-1900) along circular profile around central aster outside actin cluster vs Angle (deg) (magenta curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4d(image, at, "fascin")
