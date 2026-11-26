"""Calibrated figure reader wrapper for figS6b_data."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figS6b

def read(image, at):
    return _figS6b.read_series(image, 'data', at)
