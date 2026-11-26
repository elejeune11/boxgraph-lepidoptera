"""Calibrated figure reader wrapper for fig2f_protrusions."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig2efg

def read(image, at):
    return _fig2efg.read_panel(image, '2f')
