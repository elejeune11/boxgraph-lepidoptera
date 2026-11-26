"""Calibrated reader for Fig. 4c: Normalized fascin fluorescence intensity along dashed line across central-aster GUV vs Length (um) (magenta curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4c(image, at, "fascin")
