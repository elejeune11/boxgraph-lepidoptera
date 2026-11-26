"""Calibrated figure reader wrapper for fig1e_fascin_kinetics_fit."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1e

def read(image, at):
    return _fig1e.read_series(image, 'fit', at)
