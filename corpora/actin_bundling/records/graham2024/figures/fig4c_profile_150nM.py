"""Calibrated reader for Fig. 4C bottom line profile (150 nM Arp2/3)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig4c

def read(image, at):
    return _fig4c.read_profile(image, at, "150nM")
