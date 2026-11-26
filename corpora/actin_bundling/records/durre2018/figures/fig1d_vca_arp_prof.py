"""Calibrated reader for Fig. 1d: Normalized pyrene fluorescence I_norm (a.u.) vs time for 3 uM actin + 0.3 uM VCA + 0.3 uM Arp2/3 + 13.5 uM profilin (medium grey circles)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _curves_f1_s1_s7_s8

def read(image, at):
    return _curves_f1_s1_s7_s8.read_fig1d(image, at, "vca_arp_prof")
