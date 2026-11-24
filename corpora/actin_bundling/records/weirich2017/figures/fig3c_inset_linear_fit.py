"""Calibrated reader for Fig. 3C inset dashed linear fit line."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fig3c_inset import read_series


def read(image, at):
    return read_series(image, at, "linear_fit")
