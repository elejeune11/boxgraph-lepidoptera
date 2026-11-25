"""Calibrated reader for Fig. 3C (aspect ratio across 4 Arp2/3 & 6 uM VCA conditions)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _boxplot

def read(image, at):
    # Boxes lie below y = 3.5; brackets sit above y = 4.2
    return _boxplot.read_boxplot(image, [6, 5, 4, 3, 2, 1], [1, 2, 3, 4], "Fig3C", max_y_for_boxes=3.5)
