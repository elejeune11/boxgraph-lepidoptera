"""Calibrated heatmap slice reader for Extended Data Fig. 7e: D = 1.0 um, gamma_I = 10.0 pN."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _edfig7 import read_ed7_slice


def read(image, at):
    return read_ed7_slice(image, at, d_um=1.0, label="7e_gI=10.0pN")
