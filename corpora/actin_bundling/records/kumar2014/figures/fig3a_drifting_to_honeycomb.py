"""Calibrated reader for Fig. 3(a) 2D ARD drifting-square to oscillating-honeycomb phase boundary beta(P)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3a import read_drifting_to_honeycomb


def read(image, at):
    return read_drifting_to_honeycomb(image, at)
