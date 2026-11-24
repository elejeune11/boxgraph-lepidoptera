"""Calibrated reader for Fig. 9(a) upper blue analytical line Eq. (24) for large systems."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig9


def read(image, at):
    return _fig9.read_fig9a(image, at, "eq24_large")
