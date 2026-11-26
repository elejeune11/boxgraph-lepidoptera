"""Calibrated reader for Fig. s7: FRAP single-exponential fit (dotted red line)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _figs7

def read(image, at):
    return _figs7.read_fit(image, at)
