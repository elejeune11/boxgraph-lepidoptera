"""Calibrated reader for Main Fig. 4h: overlaid experimental AR < 1.1 power-law trendline (black solid curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4h import read_exp_fit


def read(image, at):
    return read_exp_fit(image, at)
