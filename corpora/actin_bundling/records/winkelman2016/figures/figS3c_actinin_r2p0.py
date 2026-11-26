"""Calibrated figure reader wrapper for figS3c_actinin_r2p0."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figS3c

def read(image, at):
    return _figS3c.read_series(image, 'actinin_r2p0', at)
