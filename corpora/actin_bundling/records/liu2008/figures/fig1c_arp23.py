"""Calibrated reader for Fig. 1c: Arp2/3 / Actin normalized line-scan (black curve, n = 4)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1c

def read(image, at):
    return _fig1c.read_series(image, at, "arp23")
