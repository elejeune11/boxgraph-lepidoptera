"""Calibrated reader for Fig. S4C (p. 29): normalized DIC line scan (post-phase-separation addition)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS4c(image, at, "dic")
