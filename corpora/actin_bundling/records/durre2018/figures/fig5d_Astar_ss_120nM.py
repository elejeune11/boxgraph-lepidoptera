"""Calibrated reader for Fig. 5d: Horizontal dashed steady-state guide line Arp2/3*_ss for 120 nM CP in Fig. 5d."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig5cde_s11

def read(image, at):
    return _fig5cde_s11.read_fig5d_ss(image, at, 120)
