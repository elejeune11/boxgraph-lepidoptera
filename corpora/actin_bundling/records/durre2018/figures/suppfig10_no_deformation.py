"""Calibrated reader for Supp Fig. 10: No membrane deformation vesicle percentage (%) with 5 uM NMM II at c_CP = 60 and 120 nM CP (yellow top bar segment)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_s2_s4_s10

def read(image, at):
    return _bars_s2_s4_s10.read_suppfig10(image, at, "no_deformation")
