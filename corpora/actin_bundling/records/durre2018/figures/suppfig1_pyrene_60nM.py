"""Calibrated reader for Supp Fig. 1: Pyrene fluorescence I (a.u.) vs time (0-60 min) for 3 uM actin + 60 nM CP + 13.5 uM profilin (inhibited bulk polymerization)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _curves_f1_s1_s7_s8

def read(image, at):
    return _curves_f1_s1_s7_s8.read_suppfig1(image, at)
