"""Calibrated reader for Supplemental Fig. S3 (150 nM Arp2/3): shell, ring, disc abundance %."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs3

def read(image, at):
    return _figs3.read_group(image, 2, "150nM_Arp2/3")
