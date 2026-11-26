"""Calibrated reader for Fig. 3B (p. 4): pyrene fluorescence for 1.5 uM actin alone."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3b(image, at, "actin_alone")
