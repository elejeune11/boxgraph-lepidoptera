"""Calibrated reader for Fig. 1D (p. 2): normalized DIC line scan."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig1d(image, at, "dic")
