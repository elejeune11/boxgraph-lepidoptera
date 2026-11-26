"""Calibrated reader for Fig. 2d: % Ring formation with fascin (0.1 M/M) vs GUV diameter category (1=7-12, 2=12-16, 3=>16 um) at alpha-actinin/actin = 0.1 (blue bars + SEM error bars)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig2d(image, "0.1")
