"""Calibrated figure reader wrapper for fig3h_coop_model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3h

def read(image, at):
    return _fig3h.read_series(image, 'coop', at)
