"""Calibrated reader for Fig. 3C (p. 5): V-A-CP (VASP + 2 uM actin/LatB + capping protein, blue curve)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3c_series(image, at, 'v_a_cp')
