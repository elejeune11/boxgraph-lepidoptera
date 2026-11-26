"""Calibrated reader for Supp Fig. 8: Pyrene fluorescence intensity (A.U.) vs time (min) for branched actin network polymerization at +60 nM CP."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _curves_f1_s1_s7_s8

def read(image, at):
    return _curves_f1_s1_s7_s8.read_suppfig8(image, at, 60)
