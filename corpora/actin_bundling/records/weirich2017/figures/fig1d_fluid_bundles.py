"""Calibrated reader for Fig. 1D fluid bundles phase points."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig1d import read_class


def read(image, at):
    return read_class(image, "fluid_bundles")
