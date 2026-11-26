"""Calibrated figure reader wrapper for fig3i_fascin_domains."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig3i

def read(image, at):
    return _fig3i.read_series(image, 'fascin_domains')
