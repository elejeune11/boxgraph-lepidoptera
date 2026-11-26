"""Calibrated reader for Fig. 1c: normalized Actin line-scan (red curve, n = 7)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig1c

def read(image, at):
    return _fig1c.read_series(image, at, "actin")
