"""Calibrated reader for Fig. 5E ([Arp2/3]=25 nM): filament length cumulative density function."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5_curves

def read(image, at):
    return _fig5_curves.read_fig5e(image, at, 25)
