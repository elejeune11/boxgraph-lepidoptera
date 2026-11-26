"""Calibrated reader for Fig. 3D (p. 4): excimer fluorescence for + 3 uM pLK."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3d(image, at, "plus_plk")
