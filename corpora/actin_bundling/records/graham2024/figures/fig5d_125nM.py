"""Calibrated reader for Fig. 5D ([Arp2/3]=125 nM): number of actin filaments vs time."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5_curves

def read(image, at):
    return _fig5_curves.read_fig5d(image, at, 125)
