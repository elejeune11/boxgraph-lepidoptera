"""Calibrated reader for Fig. 3C (p. 4): assembly rate 1/t_1/2 for + (pLK + pRE)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3c(image, at, "coacervate")
