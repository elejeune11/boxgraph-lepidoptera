"""Calibrated figure reader wrapper for figS2f_125nM."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figS2f

def read(image, at):
    return _figS2f.read_series(image, '125nM', at)
