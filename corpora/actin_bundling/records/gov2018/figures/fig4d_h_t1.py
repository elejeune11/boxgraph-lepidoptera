"""Calibrated reader for Fig. 4(d) top: h(x) at t = 1 s (blue, left scale)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4d import read_fig4d


def read(image, at):
    return read_fig4d(image, at, 'h', 't1')
