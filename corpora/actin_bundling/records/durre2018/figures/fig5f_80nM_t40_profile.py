"""Calibrated reader for Fig. 5f: Finite-pool cone radial barbed end density profile B(r, t=40 s) at 80 nM CP along the unoccluded left-to-center front cross-section (r = -2.5 to +0.5 um), showing the double-peaked ring wave."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5f_s14

def read(image, at):
    return _fig5f_s14.read_fig5f_80nM_t40_profile(image, at)
