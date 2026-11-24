"""Calibrated reader for Fig. 4(d) middle: m(x) at t = 2 s (green)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4d import read_fig4d


def read(image, at):
    return read_fig4d(image, at, 'm', 't2')
