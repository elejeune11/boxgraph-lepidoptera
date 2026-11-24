"""Calibrated reader for Fig. 4(d) inset top: h(x) at t = 120 s (red)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4d_inset import read_fig4d_inset


def read(image, at):
    return read_fig4d_inset(image, at, 't120')
