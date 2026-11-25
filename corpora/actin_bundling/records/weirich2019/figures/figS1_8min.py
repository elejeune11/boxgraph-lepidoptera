"""Calibrated reader for Fig. S1 — myosin cluster width vs length at 8 min (blue circles)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1 import read_series


def read(image, at):
    return read_series(image, "8min")
