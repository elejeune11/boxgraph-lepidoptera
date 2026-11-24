"""Calibrated reader for Fig. 6(a) numerical linear stability N_max vs k_f/k_0 (blue step curve)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig6


def read(image, at):
    return _fig6.read_fig6a(image, at, "Nmax")
