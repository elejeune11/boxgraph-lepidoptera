"""Calibrated reader for Fig. S9 — spherical two-droplet free energy for w/gamma = 2.5 (red dashed curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS9 import read_curve


def read(image, at):
    return read_curve(image, at, "w2p5")
