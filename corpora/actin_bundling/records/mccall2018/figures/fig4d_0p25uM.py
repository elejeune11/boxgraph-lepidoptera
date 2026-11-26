"""Calibrated reader for Fig. 4D (p. 5): 0.25 uM TMR-actin intensity line scan."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4d(image, at, "0p25uM")
