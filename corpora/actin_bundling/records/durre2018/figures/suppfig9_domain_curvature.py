"""Calibrated reader for Supp Fig. 9: Median actin domain curvature magnitude (1/um, printed unsigned/positive) at 40 nM CP (grey horizontal median bar across beeswarm scatter)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _box_scatter_f4_f6_s6_s9

def read(image, at):
    return _box_scatter_f4_f6_s6_s9.read_suppfig9(image, at, "curvature")
