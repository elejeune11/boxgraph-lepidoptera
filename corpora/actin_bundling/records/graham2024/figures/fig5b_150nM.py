"""Calibrated reader for Fig. 5B ([Arp2/3]=150 nM): radial F-actin density ratio vs subvolume radius."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5_curves

def read(image, at):
    return _fig5_curves.read_fig5b(image, at, 150)
