"""Calibrated reader for Fig. 2b: photobleached spot edge proximal to tip (black circles)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig2b

def read(image, at):
    return _fig2b.read_series(image, "proximal")
