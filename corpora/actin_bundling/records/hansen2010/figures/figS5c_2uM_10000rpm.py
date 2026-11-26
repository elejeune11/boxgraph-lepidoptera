"""Calibrated reader for Fig. S5C (p. 23): 2uM monomeric actin AUC best-fit model curve at 10000 RPM."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS5c_speed(image, at, 10000, 1.32)
