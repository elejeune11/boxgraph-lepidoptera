"""Calibrated reader for Main Fig. 3f AR < 1.1 power-law trendline (dark blue dashed)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3ef import read_power_fit


def read(image, at):
    return read_power_fit(image, at, "blue", "Fig3f AR<1.1 power fit")
