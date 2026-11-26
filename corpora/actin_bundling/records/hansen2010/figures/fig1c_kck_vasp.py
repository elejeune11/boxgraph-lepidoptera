"""Calibrated reader for Fig. 1C (p. 3): + KCK-VASP (blue pyrene fluorescence curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig1c_series(image, at, 'kck_vasp')
