"""Calibrated reader for Supplemental Fig. S5d ([VCA]=5 uM): R.A.E. vs [Arp2/3]."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs5d

def read(image, at):
    return _figs5d.read_panel(image, "VCA=5uM")
