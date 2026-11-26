"""Calibrated reader for Fig. 1b: Cumulative probability of Ring + Ring/Network (top of middle grey stacked segment; 1 minus this equals Network probability) vs GUV diameter category (1=7-12, 2=12-16, 3=>16 um)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig1b(image, "ring_plus_ringnet")
