"""Calibrated reader for Fig. 5(g): turnover-normalized ordering rate S_dot vs S_0 at kT ~ 0 (red)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5 import read_fig5g_series


def read(image, at):
    return read_fig5g_series(image, "kT_zero")
