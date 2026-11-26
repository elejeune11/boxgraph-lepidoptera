"""Calibrated figure reader wrapper for figS5e_actinin_domains_data."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figS5de

def read(image, at):
    return _figS5de.read_series(image, 'S5e', 'domains', 'data', at)
