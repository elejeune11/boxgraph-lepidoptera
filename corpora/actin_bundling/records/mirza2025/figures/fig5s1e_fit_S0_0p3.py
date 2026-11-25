"""Calibrated reader for Fig. 5-S1 III(e): black dashed Eq. (12) exponential fit for S_0 = 0.3."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig5s1e import read_fit


def read(image, at):
    return read_fit(image, at, 0.3)
