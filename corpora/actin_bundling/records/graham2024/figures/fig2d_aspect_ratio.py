"""Calibrated reader for Fig. 2D (aspect ratio vs [Arp2/3], minor axis < 3 um)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _boxplot

def read(image, at):
    return _boxplot.read_boxplot(image, [6, 5, 4, 3, 2, 1, 0], [0, 50, 100, 150], "Fig2D")
