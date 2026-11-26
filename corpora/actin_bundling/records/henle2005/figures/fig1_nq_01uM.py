"""Calibrated reader for Fig. 1 main: n_q = 0.1 uM (solid curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _henle_figs import read_fig1_main


def read(image, at):
    return read_fig1_main(image, at, "0.1uM")
