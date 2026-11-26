"""Calibrated reader for Fig. S4D (p. 29): TMR-actin RAM vs droplet radius scatter (N=117)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS4d(image)
