"""Calibrated reader for Fig. 6C: Rods abundance (%) vs [Arp2/3] (mean +/- SD, n=3)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig6c

def read(image, at):
    return _fig6c.read_stacked(image, "rods")
