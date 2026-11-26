"""Calibrated reader for Fig. s6: measured protrusion 1 elongation (11 black squares)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs6

def read(image, at):
    return _figs6.read_markers(image, "squares")
