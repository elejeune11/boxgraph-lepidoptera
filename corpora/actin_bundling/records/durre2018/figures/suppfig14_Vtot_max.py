"""Calibrated reader for Supp Fig. 14: Peak total membrane VCA density V+V*+A*+V+ (1/um^2) in finite-pool cone simulations at 40, 80, and 120 nM CP read from the three auto-scaled colorbar upper limits."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5f_s14

def read(image, at):
    return _fig5f_s14.read_sfig14_all_extrema(image, at, "max")
