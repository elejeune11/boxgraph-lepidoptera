"""Calibrated reader for Fig. 6(a) critical forward rate k_f = k_f^c (Eq. 17, right vertical black line)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig6


def read(image, at):
    return _fig6.read_fig6a(image, at, "kfc")
