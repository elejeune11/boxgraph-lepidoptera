"""Calibrated reader for Fig. 4F (150 nM Arp2/3): Shell, Ring, Disc abundance %."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig4f

def read(image, at):
    return _fig4f.read_panel(image, 2, "150nM_Arp2/3")
