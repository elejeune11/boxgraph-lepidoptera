"""Calibrated reader for Supp Fig. 4: Flat domain vesicle percentage (%) vs c_CP = 40, 60, 120, 180 nM (yellow top segment in Domains panel)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bars_s2_s4_s10

def read(image, at):
    return _bars_s2_s4_s10.read_suppfig4_domains(image, at, "flat_domain")
