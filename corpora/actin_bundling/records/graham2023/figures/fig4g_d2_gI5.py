"""Calibrated heatmap slice reader for Fig. 4g: D = 2.0 um, gamma_I = 5.0 pN."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4g import read_slice


def read(image, at):
    return read_slice(image, at, block_idx=1, d_um=2.0, label="gI=5.0pN")
