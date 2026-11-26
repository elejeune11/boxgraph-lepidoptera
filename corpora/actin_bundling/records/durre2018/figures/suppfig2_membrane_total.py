"""Calibrated reader for Supp Fig. 2: Total no-bulk / membrane-bound vesicle percentage (%) at c_CP = 0 and 20 nM CP (orange + blue stacked bars)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_s2_s4_s10

def read(image, at):
    return _bars_s2_s4_s10.read_suppfig2(image, at, "membrane_total")
