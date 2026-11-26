"""Calibrated reader for Fig. 4c: Individual domain curvature (1/um) vs domain diameter (um) scatter markers at 120 nM CP (dark grey open triangles)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _box_scatter_f4_f6_s6_s9

def read(image, at):
    return _box_scatter_f4_f6_s6_s9.read_fig4c(image, at, 120)
