"""Calibrated reader for Fig. 6c: Boxplot medians (red lines) and non-outlier whisker ranges of lamellipodium thickness (um) across 1=CP Dox (N=21), 2=no CP Dox (N=10), 3=CP no Dox (N=21)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _box_scatter_f4_f6_s6_s9

def read(image, at):
    return _box_scatter_f4_f6_s6_s9.read_fig6c(image, at)
