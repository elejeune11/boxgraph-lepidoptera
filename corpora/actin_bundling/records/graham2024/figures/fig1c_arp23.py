"""Calibrated reader for Fig. 1C: VASP + 150 nM Arp2/3 (pink) diameter histogram."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1c

def read(image, at):
    return _fig1c.read_series(image, "arp23")
