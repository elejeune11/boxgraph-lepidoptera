"""Calibrated figure reader wrapper for fig2a_fascin_pellet."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig2a

def read(image, at):
    return _fig2a.read_series(image, 'fascin')
