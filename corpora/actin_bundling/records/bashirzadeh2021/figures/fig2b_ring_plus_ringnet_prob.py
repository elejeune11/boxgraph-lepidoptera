"""Calibrated reader for Fig. 2b: Cumulative probability of Ring + Ring/Network with fascin (top of middle grey segment; 1 minus this equals Network) vs GUV diameter category (1=7-12 [printed 7-11], 2=12-16, 3=>16 um)."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def read(image, at):
    return C.read_fig2b(image, "ring_plus_ringnet")
