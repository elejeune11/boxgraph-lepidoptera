"""Calibrated reader for Supplementary Fig. 8b: Normalized actin fluorescence intensity along circular profile around star-like aster outside central aggregate vs Angle (deg) at alpha-actinin/actin=0.3, fascin/actin=0.3 (green curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig8b(image, at, "actin")
