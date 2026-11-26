"""Calibrated reader for Fig. S2C bottom (p. 26): normalized BSA line scan (BSA alone)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS2c(image, at, "bottom", "bsa")
