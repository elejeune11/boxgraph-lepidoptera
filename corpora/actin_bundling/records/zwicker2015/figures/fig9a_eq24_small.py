"""Calibrated reader for Fig. 9(a) lower blue analytical line Eq. (24) for small systems (legend prints Eq. 21)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig9


def read(image, at):
    return _fig9.read_fig9a(image, at, "eq24_small")
