"""Calibrated reader for Fig. 3E (aspect ratio across 4 Arp2/3, 750 nM VCA & CK666 conditions)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _boxplot

def read(image, at):
    return _boxplot.read_boxplot(image, [6, 5, 4, 3, 2, 1, 0], [1, 2, 3, 4], "Fig3E", max_y_for_boxes=3.5)
