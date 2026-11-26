"""Calibrated reader for Fig. 5c: Modelled barbed end density (1/um^2) vs time (0-40 s) at infinite protein pools for 80 nM CP (medium grey curve)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5cde_s11

def read(image, at):
    return _fig5cde_s11.read_fig5c(image, at, 80)
