"""Calibrated figure reader wrapper for fig5i_actinin_data."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5i

def read(image, at):
    return _fig5i.read_series(image, 'actinin', 'data', at)
