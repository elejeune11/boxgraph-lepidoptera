"""Calibrated reader for Fig. 2b: filopodium-like protrusion tip trajectory (blue circles)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig2b

def read(image, at):
    return _fig2b.read_series(image, "tip")
