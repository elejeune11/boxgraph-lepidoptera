"""Calibrated reader for Fig. 1C: VASP droplets (dark purple) diameter histogram."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1c

def read(image, at):
    return _fig1c.read_series(image, "vasp")
