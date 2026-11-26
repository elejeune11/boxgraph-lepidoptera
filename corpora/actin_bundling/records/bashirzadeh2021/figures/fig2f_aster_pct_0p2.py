"""Calibrated reader for Fig. 2f: % Aster formation (peripheral or central) with fascin (0.1 M/M) vs GUV diameter category (1=7-16, 2=>16 um) at alpha-actinin/actin = 0.2 (magenta bars + one-sided upper whisker; caption states SEM, visually ~+1 SD)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig2f(image, "0.2")
