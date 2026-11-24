"""Calibrated reader for Main Fig. 3e Discs power-law trendline (dark red dashed)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3ef import read_power_fit


def read(image, at):
    return read_power_fit(image, at, "red", "Fig3e Discs power fit")
