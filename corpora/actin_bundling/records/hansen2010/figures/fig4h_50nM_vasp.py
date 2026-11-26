"""Calibrated reader for Fig. 4H (p. 7): + 50 nM VASP (light gray outer bars) with SD error bars."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4h_series(image, '50nM_vasp')
