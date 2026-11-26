"""Calibrated reader for Supplementary Fig. 9c: Actin bundle persistence length L_p (um) at alpha-actinin/actin = 0.3 vs fascin/actin ratio (0.1 [n=26], 0.3 [n=16]; grey bars + SEM error bars)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_sfig9c(image)
