"""Calibrated reader for SI Fig. S1 (p. 6): red dashed power-law scaling line R* ~ n^(5/3)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1 import read_figS1_powerlaw


def read(image, at):
    return read_figS1_powerlaw(image, at)
