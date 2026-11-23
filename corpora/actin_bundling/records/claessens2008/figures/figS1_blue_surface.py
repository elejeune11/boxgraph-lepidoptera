"""Calibrated reader for SI Fig. S1 (p. 6): theoretical full-saturation curve (blue, surface occupation included)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figS1 import read_figS1_curve


def read(image, at):
    return read_figS1_curve(image, at, "blue_surface")
