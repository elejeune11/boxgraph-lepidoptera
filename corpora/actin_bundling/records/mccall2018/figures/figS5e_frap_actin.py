"""Calibrated reader for Fig. S5E (p. 30): normalized FRAP recovery of TMR-actin."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS5e(image, at, "actin")
