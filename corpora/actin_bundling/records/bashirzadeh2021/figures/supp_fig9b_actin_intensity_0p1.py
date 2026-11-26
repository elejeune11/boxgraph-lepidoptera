"""Calibrated reader for Supplementary Fig. 9b: Actin fluorescence intensity (a.u.) line scan across central-aster GUV vs Length (um) at fascin/actin = 0.1, alpha-actinin/actin = 0.3 (red curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig9b(image, at, "0.1")
