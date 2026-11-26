"""Calibrated reader for Fig. 2 inset: frustrated attraction, n_q = 0.1 uM (solid curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _henle_figs import read_fig2_inset


def read(image, at):
    return read_fig2_inset(image, at, "0.1uM")
