"""Calibrated reader for Supp Fig. 4: Bulk polymerization vesicle percentage (%) vs c_CP = 40, 60, 120, 180 nM (white outlined bars in Bulk Polymerization panel)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_s2_s4_s10

def read(image, at):
    return _bars_s2_s4_s10.read_suppfig4_bulk(image, at)
