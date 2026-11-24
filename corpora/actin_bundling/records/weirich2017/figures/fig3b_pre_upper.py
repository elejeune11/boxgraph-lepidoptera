"""Calibrated reader for Fig. 3B pre-coalescence upper tactoid length."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3b import read_series


def read(image, at):
    return read_series(image, at, "pre_upper")
