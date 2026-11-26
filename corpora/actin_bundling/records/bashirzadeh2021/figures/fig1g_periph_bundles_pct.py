"""Calibrated reader for Fig. 1g: Cumulative % of peripheral asters in large GUVs (>16 um) having all peripheral bundles vs alpha-actinin/actin ratio (0.1, 0.2, 0.3; 100% minus this equals lumenal bundles)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig1g(image)
