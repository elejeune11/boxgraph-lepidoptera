"""Calibrated reader for Fig. S6C (p. 31): phalloidin fluorescence step-like winding trace."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS6c(image, at)
