"""Calibrated reader for Fig. S1B (p. 25): coalescing droplet aspect ratio data + dy bars."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_figS1b(image, at, "data")
