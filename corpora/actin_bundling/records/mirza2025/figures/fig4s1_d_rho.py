"""Calibrated reader for Fig. 4-S1(d) col 1: density rho/rho_0 (Reduce gamma by 200-fold)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig4s1 import read_col1_rho


def read(image, at):
    return read_col1_rho(image, at, 'd')
