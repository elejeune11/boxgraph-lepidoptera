"""Calibrated reader for Fig. 4a: Median domain diameter (um) at c_CP = 60 and 120 nM with 5 uM NMM II (yellow horizontal median bars)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _box_scatter_f4_f6_s6_s9

def read(image, at):
    return _box_scatter_f4_f6_s6_s9.read_fig4a(image, at, "myosin")
