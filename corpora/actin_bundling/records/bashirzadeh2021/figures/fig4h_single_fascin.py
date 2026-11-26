"""Calibrated reader for Fig. 4h (top): AFINES Binder cumulant U_4 vs Time (s) for fascin alone with actin (Single component, black trajectory + 1 SD error bars across 5 simulations)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig4h(image, "fascin")
