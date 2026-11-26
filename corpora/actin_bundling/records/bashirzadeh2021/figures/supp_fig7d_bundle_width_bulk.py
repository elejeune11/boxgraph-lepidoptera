"""Calibrated reader for Supplementary Fig. 7d: Mean actin bundle width Delta (um) in bulk vs fascin/actin molar ratio (0.05, 0.1, 0.5; grey bars + SD error bars, 20 bundles per condition)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig7d(image, "bulk")
