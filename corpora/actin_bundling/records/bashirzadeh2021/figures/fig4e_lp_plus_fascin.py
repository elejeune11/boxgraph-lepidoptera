"""Calibrated reader for Fig. 4e: Actin bundle persistence length L_p (um) with fascin (0.1 M/M) vs alpha-actinin/actin ratio (0.1 [n=17], 0.2 [n=14], 0.3 [n=26]; magenta bars + symmetric error bars)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4e(image, "plus_fascin")
