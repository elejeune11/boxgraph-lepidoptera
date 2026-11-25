"""Calibrated reader for Fig. 2F (aspect ratio vs [Arp2/3], 5 uM VASP + 0.67 uM actin + 4 uM VCA)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _boxplot

def read(image, at):
    return _boxplot.read_boxplot(image, [25, 20, 15, 10, 5, 0], [0, 50, 100, 150], "Fig2F")
