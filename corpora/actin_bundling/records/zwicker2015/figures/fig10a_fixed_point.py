"""Calibrated reader for Fig. 10(a) unequal-core unique stable fixed point (R_1/w, R_2/w)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fig10


def read(image, at):
    return _fig10.read_fig10a_fp(image)
