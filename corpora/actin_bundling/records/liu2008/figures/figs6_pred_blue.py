"""Calibrated reader for Fig. s6: kinetic model predicted elongation for protrusion 2 (blue curve)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs6

def read(image, at):
    return _figs6.read_curve(image, at, "blue")
