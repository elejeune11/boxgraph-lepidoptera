"""Calibrated reader for Supp Fig. 11: Modelled barbed end density (1/um^2) vs time (0-100 s) at infinite protein pools for 40 nM CP (light grey curve)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5cde_s11

def read(image, at):
    return _fig5cde_s11.read_suppfig11(image, at, 40)
