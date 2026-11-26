"""Calibrated reader for Supp Fig. 2: Bulk polymerization vesicle percentage (%) at c_CP = 0 nM (N=24) and 20 nM (N=96) (white outlined bars)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_s2_s4_s10

def read(image, at):
    return _bars_s2_s4_s10.read_suppfig2(image, at, "bulk")
