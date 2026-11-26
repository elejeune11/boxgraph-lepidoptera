"""Calibrated reader for Fig. 3d: Cumulative % Center cluster formation with alpha-actinin+fascin vs GUV diameter bin (1=16-24 um [caption 17-24 um, n=61], 2=>24 um [n=26]; white bars)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig3d(image, "center")
