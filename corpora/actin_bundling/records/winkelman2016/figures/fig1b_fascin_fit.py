"""Calibrated figure reader wrapper for fig1b_fascin_fit."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1bc

def read(image, at):
    return _fig1bc.read_series(image, '1b', 'fit', at)
