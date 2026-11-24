"""Calibrated reader for Extended Data Fig. 9c: 15% PEG droplet size distribution (yellow)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _edfig9c import read_dist


def read(image, at):
    return read_dist(image, at, "15peg")
