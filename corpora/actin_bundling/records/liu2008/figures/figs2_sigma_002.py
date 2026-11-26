"""Calibrated reader for Fig. s2: red bundling boundary at sigma = 0.02 k_B*T/nm^2."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs2

def read(image, at):
    return _figs2.read_boundary(image, at, "sigma_002")
