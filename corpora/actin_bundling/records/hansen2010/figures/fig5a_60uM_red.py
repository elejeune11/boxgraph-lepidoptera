"""Calibrated reader for Fig. 5A (p. 8): 60uM hProI AUC trace (red speed curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig5a_curve(image, at, 'red')
