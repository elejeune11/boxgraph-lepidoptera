"""Calibrated reader for Fig. 5e: Modelled capping event rate (1/(um^2*s)) vs time (0-40 s) at infinite pools for 120 nM CP (dark crimson dashed curve)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5cde_s11

def read(image, at):
    return _fig5cde_s11.read_fig5e(image, at, "cap", 120)
