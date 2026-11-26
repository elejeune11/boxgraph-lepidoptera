"""Calibrated reader for Fig. 5A (p. 8): 90uM hProI AUC trace (green speed curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig5a_curve(image, at, 'green')
