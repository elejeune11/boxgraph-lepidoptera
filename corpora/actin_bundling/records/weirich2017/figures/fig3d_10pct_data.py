"""Calibrated reader for Fig. 3D 10% filamin relaxation time vs capillary velocity."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3d import read_series


def read(image, at):
    return read_series(image, at, "10pct")
