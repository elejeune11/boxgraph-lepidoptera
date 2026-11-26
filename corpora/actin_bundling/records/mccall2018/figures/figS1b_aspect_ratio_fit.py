"""Calibrated reader for Fig. S1B (p. 25): exponential coalescence fit curve (tau = 6.33 ms)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS1b(image, at, "fit")
