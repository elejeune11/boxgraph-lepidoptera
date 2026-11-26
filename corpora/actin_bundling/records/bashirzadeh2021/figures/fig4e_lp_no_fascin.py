"""Calibrated reader for Fig. 4e: Actin bundle persistence length L_p (um) without fascin vs alpha-actinin/actin ratio (0.1 [n=22], 0.3 [n=14]; blue bars + symmetric error bars)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4e(image, "no_fascin")
