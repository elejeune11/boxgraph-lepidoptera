"""Calibrated figure reader wrapper for fig5f_espin_fascin_fit."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5f

def read(image, at):
    return _fig5f.read_series(image, 'fit', at)
