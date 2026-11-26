"""Calibrated reader for Supp Fig. 7: TIRF fluorescence intensity profile (a.u.) across one actin domain vs spatial position (um) at 120 nM CP and t = 0 min."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _curves_f1_s1_s7_s8

def read(image, at):
    return _curves_f1_s1_s7_s8.read_suppfig7(image, at, 120, 0)
