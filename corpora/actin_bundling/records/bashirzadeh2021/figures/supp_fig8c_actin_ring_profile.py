"""Calibrated reader for Supplementary Fig. 8c: Normalized actin fluorescence intensity along dashed line across encapsulated actin ring vs Length (um) (green curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig8c(image, at, "actin")
