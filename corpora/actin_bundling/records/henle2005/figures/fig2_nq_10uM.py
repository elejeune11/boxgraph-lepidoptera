"""Calibrated reader for Fig. 2 main: finite-size ions, n_q = 10 uM (dashed curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _henle_figs import read_fig2_main


def read(image, at):
    return read_fig2_main(image, at, "10uM")
