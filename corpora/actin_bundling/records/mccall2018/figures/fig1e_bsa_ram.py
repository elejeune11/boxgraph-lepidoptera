"""Calibrated reader for Fig. 1E (p. 2): 647-BSA average RAM bars + SEM."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig1e(image, "bsa")
