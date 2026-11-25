"""Calibrated reader for Supplemental Fig. S2 (CK666 vehicle DMSO control aspect ratio)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _boxplot

def read(image, at):
    return _boxplot.read_boxplot(image, [6, 5, 4, 3, 2, 1, 0], [1, 2, 3], "FigS2", max_y_for_boxes=3.5)
