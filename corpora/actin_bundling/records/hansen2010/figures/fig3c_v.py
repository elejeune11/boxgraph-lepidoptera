"""Calibrated reader for Fig. 3C (p. 5): V (VASP alone, black curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3c_series(image, at, 'v')
