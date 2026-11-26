"""Calibrated reader for Fig. 1f: % Peripheral aster formation vs GUV diameter category (1=7-16, 2=>16 um) at alpha-actinin/actin = 0.3 (olive bars + one-sided upper error whisker; caption states SEM, visually ~+1 SD to max replicate)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig1f(image, "0.3")
