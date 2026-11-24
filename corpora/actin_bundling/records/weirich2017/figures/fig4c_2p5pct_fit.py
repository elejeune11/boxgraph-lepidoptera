"""Calibrated reader for Fig. 4C 2.5% filamin continuum model fit curve."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4c import read_series


def read(image, at):
    return read_series(image, at, "2p5pct_fit")
