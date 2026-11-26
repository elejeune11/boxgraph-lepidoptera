"""Calibrated reader for Supp Fig. 6: Boxplot medians (central grey line) and non-outlier whisker ranges of actin domains per vesicle for 1=protrusion, 2=concave domain, 3=flat domain."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _box_scatter_f4_f6_s6_s9

def read(image, at):
    return _box_scatter_f4_f6_s6_s9.read_suppfig6(image, at)
