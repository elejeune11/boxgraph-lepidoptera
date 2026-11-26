"""Calibrated reader for Fig. 4d: Alpha-actinin fluorescence intensity (a.u., right y-axis 500-1000) along circular profile around central aster outside actin cluster vs Angle (deg) (blue curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4d(image, at, "alpha_actinin")
