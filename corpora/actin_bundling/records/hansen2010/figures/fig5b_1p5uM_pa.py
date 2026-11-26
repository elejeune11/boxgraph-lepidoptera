"""Calibrated reader for Fig. 5B (p. 8): 1.5 uM profilin-actin polymerization rates across 8 conditions."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig5b_series(image, 1)
