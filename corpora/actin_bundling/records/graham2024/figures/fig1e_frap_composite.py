"""Calibrated reader for Fig. 1E: composite VASP + Arp2/3 + VCA FRAP recovery curve (green)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1e

def read(image, at):
    return _fig1e.read_curve(image, at, "composite")
